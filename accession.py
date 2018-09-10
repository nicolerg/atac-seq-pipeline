import json
import tempfile
import base64
import pdb
import encode_utils as eu
from itertools import chain
from encode_utils.connection import Connection
from google.cloud import storage


STEP_VERSION_ALIASES = {
    'default': {
        "trim-align-filter-step-v-1": '',
        "peak-call-step-v-1": '',
        "filtered-peaks-conversion-step-v-1": '',
        "signals-step-v-1": '',
        "idr-step-v-1": '',
        "overlap-peaks-step-v-1": '',
        "replicated-peaks-conversion-step-v-1": '',
        "idr-peaks-conversion-step-v-1": ''
    }
}

COMMON_METADATA = {
    'lab': '/labs/anshul-kundaje/',
    'award': 'U41HG007000'
}

RAW_FILE = {
    "aliases": [""],
    "dataset": "",
    "file_format": "",
    "flowcell_details": {
      "barcode": "",
      "flowcell": "",
      "lane": "",
      "machine": ""
    },
    "output": "",
    "paired_end": "",
    "platform": "",
    "read_length": 0,
    "replicate": "",
    "submitted_file_name": ""
}

ASSEMBLIES = ['GRCh38']


class GCBackend():
    """docstring for GCBackend"""
    def __init__(self, bucket):
        self.client = storage.Client()
        self.bucket = self.client.get_bucket(bucket)
        self.local_mapping = {}

    def blob_from_filename(self, filename):
        blob = storage.blob.Blob(self.file_path(filename), self.bucket)
        blob.reload()
        return blob

    # Returns md5sum of the file in hex
    def md5sum(self, file):
        blob = self.blob_from_filename(file)
        return self.md5_from_blob(blob)

    def size(self, file):
        blob = self.blob_from_filename(file)
        return blob.size

    # Converts base64 hash to hex format
    def md5_from_blob(self, blob):
        return base64.b64decode(blob.md5_hash).hex()

    # File path without bucket name
    def file_path(self, file):
        file_path = file.split('gs://{}/'.format(self.bucket.name))[1]
        return file_path

    # Downloads file as string
    def read_file(self, file):
        blob = self.blob_from_filename(file)
        return blob.download_as_string()

    # Downloads file to local filesystem
    def download(self, file):
        blob = self.blob_from_filename(file)
        temp_file = tempfile.NamedTemporaryFile(delete=False)
        with open(temp_file.name, 'wb'):
            blob.download_to_filename(temp_file.name)
        self.local_mapping[file] = [temp_file.name, self.md5_from_blob(blob)]
        return self.local_mapping[file]


class Analysis(object):
    """docstring for Analysis"""
    def __init__(self, metadata_json, bucket_name):
        self.files = []
        with open(metadata_json) as json_file:
            self.metadata = json.load(json_file)
        self.backend = GCBackend(bucket_name)
        self.tasks = self.make_tasks()

    # Makes instances of Task
    def make_tasks(self):
        tasks = []
        for key, value in self.metadata['calls'].items():
            for task in value:
                tasks.append(self.make_task(key, task))
        for task in tasks:
            task.output_files = self.get_or_make_files(task.outputs, task)
        # Making input files after making output files avoids creating
        # a duplicate file
        for task in tasks:
            task.input_files = self.get_or_make_files(task.inputs,
                                                      used_by_tasks=task)
        return tasks

    # Makes an instance of task with input and output GSFile instances
    def make_task(self, task_name, task):
        new_task = Task(task_name.split('.')[1], task, self)
        return new_task

    # Makes instances of GSFile from input or output section of task
    # When task=None, file is not associated with a task
    def get_or_make_files(self, section, task=None, used_by_tasks=None):
        files = []
        for key, value in section.items():
            for filename in self.extract_files(value):
                files.append(self.get_or_make_file(key,
                                                   filename,
                                                   task,
                                                   used_by_tasks))
        return files

    # Returns a GSFile object, makes a new one if one doesn't exist
    def get_or_make_file(self, key, filename, task=None, used_by_tasks=None):
        for file in self.files:
            if filename == file.filename:
                if key not in file.filekeys:
                    file.filekeys.append(key)
                if used_by_tasks not in file.used_by_tasks:
                    file.used_by_tasks.append(used_by_tasks)
                return file
        md5sum = self.backend.md5sum(filename)
        size = self.backend.size(filename)
        new_file = GSFile(key, filename, md5sum, size, task, used_by_tasks)
        self.files.append(new_file)
        return new_file

    # Cromwell workflow id
    @property
    def workflow_id(self):
        return self.metadata['labels']['cromwell-workflow-id']

    # Files in the 'outputs' of the metadata that are
    # used for filtering out intermediate outputs
    @property
    def outputs_whitelist(self):
        return list(self.extract_files(self.metadata['outputs']))

    # Files in the 'inputs' of the metadata that are
    # used for filtering out intermediate inputs
    @property
    def inputs_whitelist(self):
        return list(self.extract_files(self.metadata['inputs']))

    # Extracts file names from dict values
    def extract_files(self, outputs):
        if (isinstance(outputs, str)
                and 'gs://' in outputs
                and self.backend.bucket.name in outputs):
            yield outputs
        elif isinstance(outputs, list):
            for item in outputs:
                yield from self.extract_files(item)
        elif isinstance(outputs, dict):
            for key, values in outputs.items():
                yield from self.extract_files(values)

    def get_tasks(self, task_name):
        tasks = []
        for task in self.tasks:
            if task_name == task.task_name:
                tasks.append(task)
        return tasks

    def get_file(self, filekey=None, filename=None):
        if filekey:
            for file in self.files:
                if filekey in file.filekeys:
                    return file
        if filename:
            for file in self.files:
                if filename == file.filename:
                    return file

    @property
    def raw_fastqs(self):
        fastqs = []
        for file in self.files:
            if 'fastqs' in file.filekeys and file.task is None:
                fastqs.append(file)
        return fastqs


class Task(object):
    """docstring for Task"""
    def __init__(self, task_name, task, analysis):
        super().__init__()
        self.task_name = task_name
        self.input_files = []
        self.output_files = []
        self.inputs = task['inputs']
        self.outputs = task['outputs']
        self.docker_image = task.get('dockerImageUsed', None)
        self.analysis = analysis


class GSFile(object):
    """docstring for File"""
    def __init__(self, key, name, md5sum, size, task=None, used_by_tasks=None):
        super().__init__()
        self.filename = name
        self.filekeys = [key]
        self.task = task
        self.used_by_tasks = [used_by_tasks]
        self.md5sum = md5sum
        self.size = size


class Accession(object):
    """docstring for Accession"""
    def __init__(self, accession_id, metadata_json, bucket_name, server):
        super(Accession, self).__init__()
        self.accession_id = accession_id
        self.analysis = Analysis(metadata_json, bucket_name)
        self.backend = self.analysis.backend
        self.conn = Connection(server)

    def accession_fastqs(self):
        pass

    def file_at_portal(self, file):
        md5sum = self.backend.md5sum(file)
        search_param = [('md5sum', md5sum), ('type', 'File')]
        encode_file = self.conn.search(search_param)
        if len(encode_file) > 0:
            return encode_file[0]
        return False

    def raw_fastq_inputs(self, file):
        if not file.task and 'fastqs' in file.filekeys:
            yield file
        if file.task:
            for input_file in file.task.input_files:
                yield from self.raw_fastq_inputs(input_file)

    def raw_files_accessioned(self):
        for file in self.analysis.raw_fastqs:
            if not self.file_at_portal(file.filename):
                return False
        return True


    def accession_file(self, encode_file, gs_file):
        #local_file = self.backend.download(gs_file.filename)
        pdb.set_trace()
        encode_posted_file = self.conn.post(encode_file)
        return encode_posted_file


    def get_or_make_step_run(self, lab_prefix, run_name, step_version, task_name):
        docker_tag = self.analysis.get_tasks(task_name)[0].docker_image.split(':')[1]
        payload = {'aliases': ["{}:{}-{}".format(lab_prefix, run_name, docker_tag)],
                   'status': 'released',
                   'analysis_step_version': step_version}
        payload[Connection.PROFILE_KEY] = 'analysis_step_runs'
        return self.conn.post(payload)

    def accession_qc_objects(self):
        pass

    def make_qc_object(self):
        pass

    @property
    def assembly(self):
        assembly = [reference
                    for reference
                    in ASSEMBLIES
                    if reference
                    in self.analysis.get_tasks('read_genome_tsv')[0].outputs.get(
                        'genome', {}).get('ref_fa', '')]
        return assembly[0] if len(assembly) > 0 else ''

    def make_alignment_bam(self, file, step_run):
        encode_fastqs = [self.file_at_portal(fastq.filename)
                         for fastq
                         in set(list(self.raw_fastq_inputs(file)))]
        derived_from = ['/files/{}/'.format(obj['accession'])
                        for obj
                        in encode_fastqs]
        technical_replicates = list(set(chain(*[obj['technical_replicates']
                                        for obj
                                        in encode_fastqs])))
        biological_replicates = list(set(chain(*[obj['biological_replicates']
                                         for obj
                                         in encode_fastqs])))
        filename_for_alias = file.filename.split('gs://')[-1]
        alignment_bam = {
            'file_format':              'bam',
            'output_type':              'alignments',
            'file_type':                'bam',
            'output_category':          'alignment',
            'status':                   'released',
            'aliases':                  ['anshul-kundaje:{}'.format(filename_for_alias)],
            'assembly':                 self.assembly,
            'submitted_file_name':      file.filename.split('gs://')[-1],
            'biological_replicates':    biological_replicates,
            'technical_replicates':     technical_replicates,
            'dataset':                  encode_fastqs[0].get('dataset'),
            'step_run':                 step_run.get('@id'),
            'derived_from':             derived_from,
            'file_size':                file.size,
            'md5sum':                   file.md5sum}
        alignment_bam[Connection.PROFILE_KEY] = 'file'
        alignment_bam.update(COMMON_METADATA)
        return alignment_bam

    def accession_alignment_outputs(self,
                                    task_name='filter',
                                    filekey='nodup_bam'):
        alignment_bams = []
        tasks = self.analysis.get_tasks(task_name)
        for task in tasks:
            for bam in [file
                        for file
                        in task.output_files
                        if filekey in file.filekeys]:
                print(task_name)
                step_run = self.get_or_make_step_run(
                    'anshul-kundaje',
                    'atac-seq-trim-align-filter-step-run-single-rep-v1',
                    'anshul-kundaje:atac-seq-trim-align-filter-step-version-single-rep-v1',
                    task_name)
                encode_bam_file = self.accession_file(self.make_alignment_bam(
                    bam, step_run), bam)
                alignment_bams.append(encode_bam_file)
        return alignment_bams

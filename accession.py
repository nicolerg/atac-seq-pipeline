import json
import tempfile
import pdb
import operator
import encode_utils as eu
from itertools import chain
from functools import reduce
from base64 import b64encode, b64decode
from encode_utils.connection import Connection
from google.cloud import storage


COMMON_METADATA = {
    'lab': '/labs/anshul-kundaje/',
    'award': 'U41HG007000'
}


ASSEMBLIES = ['GRCh38', 'mm10']


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
        return b64decode(blob.md5_hash).hex()

    # File path without bucket name
    def file_path(self, file):
        file_path = file.split('gs://{}/'.format(self.bucket.name))[1]
        return file_path

    # Downloads file as string
    def read_file(self, file):
        blob = self.blob_from_filename(file)
        return blob.download_as_string()

    # Read json file
    def read_json(self, file):
        return json.loads(self.read_file(file.filename).decode())

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
                if used_by_tasks and used_by_tasks not in file.used_by_tasks:
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

    def get_files(self, filekey=None, filename=None):
        files = []
        if filekey:
            for file in self.files:
                if filekey in file.filekeys:
                    files.append(file)
        if filename:
            for file in self.files:
                if filename == file.filename:
                    files.append(file)
        return list(set(files))

    @property
    def raw_fastqs(self):
        fastqs = []
        for file in self.files:
            if 'fastqs' in file.filekeys and file.task is None:
                fastqs.append(file)
        return fastqs

    # Search the Analysis hirearchy up for a file matching filekey
    # Returns generator object, access with next() or list()
    def search_up(self, task, task_name, filekey, inputs=False):
        if task_name == task.task_name:
            if inputs:
                for file in task.input_files:
                    if filekey in file.filekeys:
                        yield file
            else:
                for file in task.output_files:
                    if filekey in file.filekeys:
                        yield file
        for task_item in set(map(lambda x: x.task, task.input_files)):
            if task_item:
                yield from self.search_up(task_item, task_name, filekey, inputs)

    # Search the Analysis hirearchy down for a file matching filekey
    # Returns generator object, access with next()
    def search_down(self, task, task_name, filekey):
        if task_name == task.task_name:
            for file in task.output_files:
                if filekey in file.filekeys:
                    yield file
        for task_item in set(reduce(operator.concat,
                                    map(lambda x: x.used_by_tasks,
                                        task.output_files))):
            if task_item:
                yield from self.search_down(task_item, task_name, filekey)


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
        self.used_by_tasks = [used_by_tasks] if used_by_tasks else []
        self.md5sum = md5sum
        self.size = size

    # Depends on all other tasks and files having finished initializing
    # Returns lisf of files
    def derived_from(self, filekey=None):
        if not filekey:
            return self.task.input_files
        else:
            return list(filter(lambda x: filekey in x.filekeys,
                               self.task.input_files))


class Accession(object):
    """docstring for Accession"""
    def __init__(self, accession_id, metadata_json, bucket_name, server):
        super(Accession, self).__init__()
        self.accession_id = accession_id
        self.analysis = Analysis(metadata_json, bucket_name)
        self.backend = self.analysis.backend
        self.conn = Connection(server)
        self.new_files = []

    def accession_fastqs(self):
        pass

    def wait_for_portal(self):
        pass

    def file_at_portal(self, file):
        self.wait_for_portal()
        md5sum = self.backend.md5sum(file)
        search_param = [('md5sum', md5sum), ('type', 'File')]
        encode_file = self.conn.search(search_param)
        if len(encode_file) > 0:
            return self.conn.get(encode_file[0].get('accession'))

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
        file_exists = self.file_at_portal(gs_file.filename)
        if not file_exists:
            local_file = self.backend.download(gs_file.filename)[0]
            encode_file['submitted_file_name'] = local_file
            encode_posted_file = self.conn.post(encode_file)
            # self.conn.upload_file(file_id=encode_posted_file.get('accession'),
            #                       file_path=local_file)
            self.new_files.append(encode_posted_file)
            return encode_posted_file
        return file_exists

    def patch_file(self, encode_file, new_properties):
        new_properties[self.conn.ENCID_KEY] = encode_file.get('accession')
        self.conn.patch(new_properties)

    def get_or_make_step_run(self, lab_prefix, run_name, step_version, task_name):
        docker_tag = self.analysis.get_tasks(task_name)[0].docker_image.split(':')[1]
        payload = {'aliases': ["{}:{}-{}".format(lab_prefix, run_name, docker_tag)],
                   'status': 'released',
                   'analysis_step_version': step_version}
        payload[Connection.PROFILE_KEY] = 'analysis_step_runs'
        return self.conn.post(payload)

    @property
    def assembly(self):
        assembly = [reference
                    for reference
                    in ASSEMBLIES
                    if reference
                    in self.analysis.get_tasks('read_genome_tsv')[0].outputs.get(
                        'genome', {}).get('ref_fa', '')]
        return assembly[0] if len(assembly) > 0 else ''

    @property
    def lab_pi(self):
        return COMMON_METADATA['lab'].split('/labs/')[1].split('/')[0]

    @property
    def dataset(self):
        return '/experiments/{}/'.format(self.accession_id)

    def file_from_template(self,
                           file,
                           file_format,
                           output_type,
                           step_run,
                           derived_from,
                           dataset):
        file_name = file.filename.split('gs://')[-1].replace('/', '-')
        obj = {
            'status':               'uploading',
            'aliases':              ['{}:{}'.format(self.lab_pi, file_name)],
            'file_format':          file_format,
            'output_type':          output_type,
            'assembly':             self.assembly,
            'submitted_file_name':  file.filename.split('gs://')[-1],
            'dataset':              dataset,
            'step_run':             step_run.get('@id'),
            'derived_from':         derived_from,
            'file_size':            file.size,
            'md5sum':               file.md5sum
        }
        if (file_format in ['bed', 'bigBed']
                and output_type
                in ['filtered peaks',
                    'conservative idr thresholded peaks',
                    'optimal idr thresholded peaks',
                    'conservative replicated peaks',
                    'replicated peaks']):
            obj['file_format_type'] = 'narrowPeak'
        obj[Connection.PROFILE_KEY] = 'file'
        obj.update(COMMON_METADATA)
        return obj

    # Returns list of accession ids of files on portal or recently accessioned
    def get_derived_from(self, file, task_name, filekey, inputs=False):
        derived_from_files = list(set(list(self.analysis.search_up(file.task,
                                                                   task_name,
                                                                   filekey,
                                                                   inputs))))
        encode_files = [self.file_at_portal(gs_file.filename)
                        for gs_file
                        in derived_from_files]
        accessioned_files = encode_files + self.new_files
        derived_from_accession_ids = []
        for gs_file in derived_from_files:
            for encode_file in accessioned_files:
                if gs_file.md5sum == encode_file.get('md5sum'):
                    derived_from_accession_ids.append(encode_file.get('accession'))
        if len(derived_from_accession_ids) != len(derived_from_files):
            raise Exception('Missing derived_from files on the portal')
        return ['/files/{}/'.format(accession_id)
                for accession_id in derived_from_accession_ids]

    # File object to be accessioned
    # inputs=True will search for input fastqs in derived_from
    def make_file_obj(self, file, file_format, output_type, step_run,
                      derived_from_task, derived_from_filekey, inputs=False):
        derived_from = self.get_derived_from(file,
                                             derived_from_task,
                                             derived_from_filekey,
                                             inputs)
        return self.file_from_template(file,
                                       file_format,
                                       output_type,
                                       step_run,
                                       derived_from,
                                       self.dataset)

    def attach_flagstat_qc_to(self, encode_bam_file, gs_file):
        qc = self.backend.read_json(self.analysis.get_files('qc_json')[0])
        flagstat_qc = qc['nodup_flagstat_qc'][int(encode_bam_file.get('biological_replicates')[0]) - 1]
        for key, value in flagstat_qc.items():
            if '_pct' in key:
                flagstat_qc[key] = '{}%'.format(value)
        step_run = encode_bam_file.get('step_run')
        if isinstance(step_run, str):
            step_run_id = step_run
        elif isinstance(step_run, dict):
            step_run_id = step_run.get('@id')
        flagstat_qc.update({
            'step_run':             step_run_id,
            'quality_metric_of':    [encode_bam_file.get('@id')],
            'status':               'released'})
        flagstat_qc.update(COMMON_METADATA)
        flagstat_qc[Connection.PROFILE_KEY] = 'samtools-flagstats-quality-metric'
        posted_qc = self.conn.post(flagstat_qc, require_aliases=False)
        return posted_qc

    def attach_cross_correlation_qc_to(self, encode_bam_file, gs_file):
        qc = self.backend.read_json(self.analysis.get_files('qc_json')[0])
        plot_pdf = next(self.analysis.search_down(gs_file.task,
                                                  'xcor',
                                                  'plot_pdf'))
        read_length_file = next(self.analysis.search_up(gs_file.task,
                                                        'bowtie2',
                                                        'read_len_log'))
        read_length = int(self.backend.read_file(read_length_file.filename).decode())
        xcor_qc = qc['xcor_score'][int(encode_bam_file.get('biological_replicates')[0]) - 1]
        pbc_qc = qc['pbc_qc'][int(encode_bam_file.get('biological_replicates')[0]) - 1]
        if isinstance(step_run, str):
            step_run_id = step_run
        elif isinstance(step_run, dict):
            step_run_id = step_run.get('@id')
        xcor_object = {
            'NRF':                  pbc_qc['NRF'],
            'PBC1':                 pbc_qc['PBC1'],
            'PBC2':                 pbc_qc['PBC2'],
            'NSC':                  xcor_qc['NSC'],
            'RSC':                  xcor_qc['RSC'],
            'sample size':          xcor_qc['num_reads'],
            "fragment length":      xcor_qc['est_frag_len'],
            "quality_metric_of":    [encode_bam_file.get('@id')],
            "step_run":             step_run_id,
            "paired-end":           self.analysis.metadata['inputs']['atac.paired_end'],
            "read length":          read_length,
            "status":               "released",
            "cross_correlation_plot": self.get_attachment(plot_pdf, 'application/pdf')
        }
        # "cross_correlation_plot": self.get_attachment(plot_pdf, 'application/pdf')

        xcor_object.update(COMMON_METADATA)
        xcor_object[Connection.PROFILE_KEY] = 'complexity-xcorr-quality-metrics'
        posted_qc = self.conn.post(xcor_object, require_aliases=False)
        return posted_qc

    def file_has_qc(self, bam, qc):
        for item in bam['quality_metrics']:
            if item['@type'][0] == qc['@type'][0]:
                return True
        return False

    def get_attachment(self, gs_file, mime_type):
        contents = self.backend.read_file(gs_file.filename)
        contents = b64encode(contents)
        if type(contents) is bytes:
            # The Portal treats the contents as string "b'bytes'"
            contents = str(contents).replace('b', '', 1).replace('\'', '')
        pdb.set_trace()
        obj = {
            'type': mime_type,
            'download': gs_file.filename.split('/')[-1],
            'href': 'data:{};base64,{}'.format(mime_type,
                                               contents)
        }
        return obj

    def accession_step(self, single_step_params):
        step_run = self.get_or_make_step_run(
            self.lab_pi,
            single_step_params['dcc_step_run'],
            single_step_params['dcc_step_version'],
            single_step_params['wdl_task_name'])
        accessioned_files = []
        for task in self.analysis.get_tasks(task_name):
            for file_params in single_step_params['wdl_files']:
                for wdl_file in [file
                                 for file
                                 in task.output_files
                                 if file_params['filekey']
                                 in file.filekeys]:
                    obj = self.make_file_obj(wdl_file,
                                             file_params['file_format'],
                                             file_params['output_type'],
                                             step_run,
                                             file_params['derived_from_task'],
                                             file_params['derived_from_filekey'])
                    encode_file = self.accession_file(obj, wdl_file)
                    # Need to move QC object adding after all files are accessioned
                    # if not list(filter(lambda x: 'SamtoolsFlagstatsQualityMetric'
                    #                              in x['@type'],
                    #                    encode_file['quality_metrics'])):
                    #     self.attach_flagstat_qc_to(encode_file, bam)
                    # if not list(filter(lambda x: 'ComplexityXcorrQualityMetric'
                    #                              in x['@type'],
                    #                    encode_file['quality_metrics'])):
                    #     self.attach_cross_correlation_qc_to(encode_file, bam)

                    accessioned_files.append(encode_file)
        return accessioned_files

    def accession_steps(self, steps_and_params_json):
        for step in steps_and_parameters_json:
            self.accession_step(step)


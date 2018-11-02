workflow test_work {
	Array[Object] steps

	call test { input :
		test_ = steps
	}


}

task test {
	Array[Object] test_
	File test_file

	command {
		cat ${test_file} > test.json
	}

	output {
		File json = glob('*.json')[0]
	}
}
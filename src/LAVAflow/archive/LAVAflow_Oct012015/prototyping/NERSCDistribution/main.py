#!/usr/bin/python

import os, sys, string, glob, math, time

def expand_arguments(inputDir, runInfo, argList):

	# print "\n\n"
	# print "In expand_arguments"

	# for r in runInfo:
	# 	print r

	newArgs = []

	for a in argList:
		# print "argument: %s"%a

		if a=='pltfile':
			# print "Expanding pltfile"
			newArgs = newArgs + [runInfo[0]]
		elif a=='prtfile':
			# print "Expanding pltfile"
			newArgs = newArgs + [string.replace(runInfo[0],'.plt.','.prt.')]
		elif a=='expfile':
			# print "Expanding expfile"
			newArgs = newArgs + ["%s%s"%(runInfo[1],'.expstat')]
		elif a=='scrfile':
			# print "Expanding scrfile"
			newArgs = newArgs + [os.path.join(inputDir,'scr_%s'%runInfo[2])]
		elif a=='outputbase':
			# print "Expanding outputbase"
			newArgs = newArgs + [runInfo[1]]
		elif a=='modelid':
			# print "Expanding modelid"
			newArgs = newArgs + [runInfo[2]]
		elif a=='modeldim':
			# print "Expanding modeldim"
			newArgs = newArgs + [runInfo[3]]
		else:
			print "ERROR: Unknown executable argument \"%s\""%a
			sys.exit(-1);


	# print newArgs

	# Return new argument list
	return newArgs



	# print "\n\n"


def parse_model_list(fileName):

	lines = open(fileName,'r').readlines(1000)
	modelID = []
	modelDim = []

	for l in lines:

		lineTrimmed = l.strip()

		# Make sure the length of the line isn't zero
		if(lineTrimmed.__len__()==0):
			continue

		# Make sure we're not a comment
		if(lineTrimmed[0]=='#'):
			continue


		# Split the line
		lineSplit = string.split(lineTrimmed)

		# Put the modelID and the model dimension into the appropriate lists
		modelID = modelID + [lineSplit[0]]
		modelDim = modelDim + [lineSplit[1]]

	# zip them into a dictionary
	res = dict(zip(modelID,modelDim))

	# return the dictionary of modelID->modelDim
	return res

def get_files_to_process(modelDict, inputDirFormat, outputDirFormat, analysisFileFormats, dataFileGlobIdentifier):
	# For every model to analyze, get every file containing dataFileIdentifiers
	inputFiles = []
	for modelID in modelDict.keys():
		modelDim = modelDict[modelID]

		# Replace the identifiers in the input directory format to obtain the current input directory
		inputDir = string.replace(inputDirFormat,'MODELID',modelID)
		inputDir = string.replace(inputDir,'MODELDIM',modelDim)
		# Replace the identifiers in the input directory format to obtain the current input directory
		outputDir = string.replace(outputDirFormat,'MODELID',modelID)
		outputDir = string.replace(outputDir,'MODELDIM',modelDim)

		print 'Input directory = %s'%inputDir

		if(not os.path.isdir(inputDir)):
			print "Input Directory does not exist!"
		if(not os.path.isdir(outputDir)):
			print "Output Directory does not exist!"
		
		# Get all possible input files
		inputFileFormat = os.path.join(inputDir,dataFileGlobIdentifier)
		inputFilesTmp = glob.glob(inputFileFormat)
		inputFilesTmp.sort()

		# Determine if the input files have already been processed based on the analysis file formats
		allOutputsExist = True
		for inFile in inputFilesTmp:

			for fileExt in analysisFileFormats:
				fileToCheck = string.replace(inFile,'.vtk','.%s'%fileExt)
				fileToCheck = string.replace(fileToCheck,inputDir,outputDir)
				
				#print fileToCheck
				
				if(not os.path.isfile(fileToCheck)):
					allOutputsExist = False
					break

			# If some or all of the output files are missing, add the input file to the list that should be processed
			if(not allOutputsExist):
				outputFileBase = string.replace(inFile,'.vtk','')
				outputFileBase = string.replace(outputFileBase,inputDir,outputDir)

				inputFiles = inputFiles + [[inFile,outputFileBase,modelID,modelDim]]
		
		if(not allOutputsExist):
			print "Not all outputs exist for directory %s"%inputDir


	# sort the list alphabetically
	#inputFiles.sort()

	# return
	return inputFiles

def chunks(a, n):
    k, m = len(a) / n, len(a) % n
    return list(a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in xrange(n))

def print_chunks(processingChunks):

	for ch in processingChunks:
		print ch

def main():
	fileName = 'modelsToAnalyze.txt'
	modelDict = parse_model_list(fileName)
	print modelDict

#	inputDirFormat = '/data1/sne/HOTB/MODELDIMd/MODELID'
#	outputDirFormat = '/data1/sne/HOTB/MODELDIMd/MODELID/output'
#	executable = '/home/tah09e/code/workspace/LAVAflow/build/Tests/ExplosionStatistics2d'

	inputDirFormat = '/global/homes/t/tah09e/project/sne/HOTB/MODELDIMd/MODELID'
	outputDirFormat = '/global/homes/t/tah09e/scratch/sneAnalysis/MODELDIMd/MODELID'

	# inputDirFormat = '/home/tah09e/code/workspace/LAVAflow/prototyping/NERSCDistribution/input/MODELDIMd/MODELID'
	# outputDirFormat = '/home/tah09e/code/workspace/LAVAflow/prototyping/NERSCDistribution/MODELDIMd/MODELID'

#	executable = '/global/homes/t/tah09e/scratch/sneAnalysis/executables/ExplosionStatistics2d'
#	executable = '/global/homes/t/tah09e/scratch/sneAnalysis/executables/ExplosionStatistics3d'
#	analysisFileFormats = ['expstat']
#	exec_args = ['pltfile','outputbase','scrfile']

#	executable = '/global/homes/t/tah09e/scratch/sneAnalysis/executables/Flux2d'
#	executable = '/global/homes/t/tah09e/scratch/sneAnalysis/executables/Flux3d'
#	analysisFileFormats = ['shell','surf']
#	exec_args = ['pltfile','outputbase']

#	executable = '/global/homes/t/tah09e/scratch/sneAnalysis/executables/ReynoldsStresses2d'
#	executable = '/global/homes/t/tah09e/scratch/sneAnalysis/executables/ReynoldsStresses3d'
#	analysisFileFormats = ['reynoldsstresses']
#	exec_args = ['pltfile','outputbase']

#	executable = '/global/homes/t/tah09e/scratch/sneAnalysis/executables/ShockSurfaceDecomposition2d'
#	executable = '/global/homes/t/tah09e/scratch/sneAnalysis/executables/ShockSurfaceDecomposition3d'
#	analysisFileFormats = ['shockdecomp']
#	exec_args = ['pltfile','outputbase']

#	executable = '/global/homes/t/tah09e/scratch/sneAnalysis/executables/RichardsonNumber2d'
#	executable = '/global/homes/t/tah09e/scratch/sneAnalysis/executables/RichardsonNumber3d'
#	analysisFileFormats = ['bruntvaisala']
#	exec_args = ['pltfile','outputbase']

	executable = '/global/homes/t/tah09e/scratch/sneAnalysis/executables/ProjectDrafts2d'
#	executable = '/global/homes/t/tah09e/scratch/sneAnalysis/executables/ProjectDrafts3d'
	analysisFileFormats = ['cluster']
	exec_args = ['pltfile','outputbase']


#	executable = '/global/homes/t/tah09e/scratch/sneAnalysis/executables/SphericalSpectra2d'
#	executable = '/global/homes/t/tah09e/scratch/sneAnalysis/executables/SphericalSpectra3d'
#	analysisFileFormats = ['spectra']
#	exec_args = ['pltfile','outputbase','expfile']

#	executable = '/global/homes/t/tah09e/scratch/sneAnalysis/executables/EulerianStructure2d'
#	executable = '/global/homes/t/tah09e/scratch/sneAnalysis/executables/EulerianStructure3d'
#	analysisFileFormats = ['evsf']
#	exec_args = ['pltfile','prtfile','outputbase','expfile']


	dataFileGlobIdentifier = '*plt*vtk'
	nJobs = 1
	qsubDir = '/global/homes/t/tah09e/scratch/sneAnalysis/qsubFiles'
	qsubDir = '/home/tah09e/code/workspace/LAVAflow/prototyping/NERSCDistribution/qsubFiles'
	writeQsubFiles = True


	# Get all of the data required to process files with LAVAflow
	# resulting format is [inFile,outputFileBase,modelID,modelDim]
	inputFiles = get_files_to_process(modelDict, inputDirFormat, outputDirFormat, analysisFileFormats, dataFileGlobIdentifier)

	# for f in inputFiles:
	# 	print f

	# Split the files to be processed into chunks
	nFilesToProcess = len(inputFiles)	

	print 'Number of files to process: %i'%len(inputFiles)
	print 'Number of jobs: %i'%nJobs


	if(nJobs>nFilesToProcess):
		print 'Number of jobs to run is greater than the number of files to process!'
		return

	processingChunks = chunks(list(i for i in range(0,nFilesToProcess)),nJobs)

#	for c in processingChunks:
#		print c

	chunkSizes = []
	for c in processingChunks:
		chunkSizes = chunkSizes + [len(c)]

#	print 'Min/Mean/Max chunk sizes: %i %f %i'%(min(chunkSizes),float(sum(chunkSizes))/len(chunkSizes) if len(chunkSizes) > 0 else float('nan'), max(chunkSizes))

	
	if(writeQsubFiles):
		# Generate qsub files
		qsubBase = """#PBS -q serial
#PBS -l walltime=00:10:00
#PBS -l pvmem=6GB
#PBS -N analysis
#PBS -e analysis.$PBS_JOBID.err
#PBS -o analysis.$PBS_JOBID.out
#PBS -m n
#PBS -V


"""
		timestamp = int(time.mktime(time.gmtime()))
		for i in range(0,len(processingChunks)):
			qsubFileName = os.path.join(qsubDir,'qsub_%s_%03i'%(timestamp,i+1))
			print qsubFileName
			
			# open the qsub file		
			fid = open(qsubFileName,'w')
			# write the header information
			fid.write(qsubBase)
			# write the call to the executable with the proper arguments for each file to process
			for j in processingChunks[i]:
				# resulting format is [inFile,outputFileBase,modelID,modelDim]
				runInfo = inputFiles[j]
				print runInfo



				inputDir = string.replace(inputDirFormat,'MODELID',runInfo[2])
				inputDir = string.replace(inputDir,'MODELDIM',runInfo[3])
				scrFile = os.path.join(inputDir,'scr_%s'%runInfo[2])
				
				args = expand_arguments(inputDir, runInfo, exec_args)
				executableLine = '%s '%executable

				for a in args:
					executableLine = executableLine + ' ' + a

				print "Executable line: %s"%executableLine

				# if(useScr):				
				# 	executableLine = '%s %s %s %s\n'%(executable,runInfo[0],runInfo[1],scrFile)
				# else:
				# 	executableLine = '%s %s %s\n'%(executable,runInfo[0],runInfo[1])
				fid.write(executableLine)

			fid.close()








if __name__ == "__main__":
	main()

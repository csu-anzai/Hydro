
import string
import os
import glob
import subprocess

def traverseTree(rootDir, outputRootDir, executable):
	for root, dirs, files in os.walk(rootDir, topdown=True): # Walk directory tree

		files = sorted(files)


		for f in files:
			#print f
			#print root
			
			# Get the model ID to find the scr file
			splitRoot = root.split('/')
			modelID = splitRoot[-1]
			#print "Model ID: %s"%modelID
		
			srcFile = os.path.join(root,"scr_%s"%modelID)
			
			#print "src file: %s"%srcFile
			
			if( string.find(f, "plt") > -1 ):
				#print "Found a checkpoint file!"
				inputFile = os.path.join(root,f)
				outputFile = inputFile.replace(rootDir,outputRootDir)
				
				statsExists = os.path.isfile("%s.%s"%(outputFile,"expstat"))
				if( statsExists ):
					continue

				
				

				subprocess.call([executable,inputFile,outputFile,srcFile])



#dataDir = "/home/tah09e/data/sne/data/HOTB/2d/b155d2a3s123SG1MHwr1LB"
#outputDir = "/home/tah09e/data/sne/output/HOTB/2d/b155d2a3s123SG1MHwr1LB"
dataDir = "/home/tah09e/data/sne/data"
outputDir = "/home/tah09e/data/sne/output"
exe = "/home/tah09e/code/workspace/LAVAflow/prototyping/TraverseDirectories/helloWorld.exe"
exe = "/home/tah09e/code/workspace/LAVAflow/build/Tests/ExplosionStatistics2d"

print "Hello world!"
traverseTree(dataDir, outputDir, exe)

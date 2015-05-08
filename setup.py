#!/usr/bin/env python
"""Description: 
Setup script for CRISPResso -- Software pipeline for the analysis of CRISPR-Cas9 genome editing outcomes from deep sequencing data
@status:  beta
@version: $Revision$
@author:  Luca Pinello
@contact: lpinello@jimmy.harvard.edu
"""

import os
import sys
import re
from setuptools import setup, Extension
import os
import pickle as cp
import glob
import subprocess as sb
import sys
import platform 
import shutil
from os.path import expanduser
import urllib


try: 
	INSTALLATION_PATH='%s/CRISPresso_dependencies' % os.environ['HOME']
except:
	INSTALLATION_PATH='%s/CRISPresso_dependencies' % os.environ['HOME']

BIN_FOLDER=os.path.join(INSTALLATION_PATH,'bin')

def main():
	if float(sys.version[:3])<2.7 or float(sys.version[:3])>=2.8:
		sys.stderr.write("CRITICAL: Python version must be 2.7!\n")
		sys.exit(1)


	version = re.search(
    	'^__version__\s*=\s*"(.*)"',
    	open('CRISPResso/bootstrap.py').read(),
    	re.M
    	).group(1)
	
	
	if float(sys.version[:3])<2.6 or float(sys.version[:3])>=2.8:
    		sys.stderr.write("ERROR: Python version must be 2.6 or 2.7!\n")
    		sys.exit(1)

	setup(
		  version=version,
          name = "CRISPResso",
    	  packages = ["CRISPResso"],
    	  package_dir={'CRISPResso': 'CRISPResso'},
          package_data={'CRISPResso': ['data/*']},     
    	  entry_points = {
        	"console_scripts": ['CRISPResso = CRISPResso.bootstrap:main']
           },
          description="Software pipeline for the analysis of CRISPR-Cas9 genome editing outcomes from deep sequencing data",
          author='Luca Pinello',
          author_email='lpinello@jimmy.harvard.edu',
          url='http://github.com/lucapinello/CRISPResso',
  
          classifiers=[
              'Development Status :: 4 - Beta',
              'Environment :: Console',
              'Intended Audience :: Developers',
              'Intended Audience :: Science/Research',              
              'License :: OSI Approved :: BSD License',
              'Operating System :: MacOS :: MacOS X',
              'Operating System :: POSIX',
              'Topic :: Scientific/Engineering :: Bio-Informatics',
              'Programming Language :: Python',
              ],
          install_requires=[
              'numpy>=1.6',
              'pandas>=0.15',
              'matplotlib>=1.3',
              'biopython>=1.6.5'
              ],

          )




def query_yes_no(question, default="yes"):
	valid = {"yes":True,   "y":True,  "ye":True,
			 "no":False,     "n":False}
	if default == None:
		prompt = " [y/n] "
	elif default == "yes":
		prompt = " [Y/n] "
	elif default == "no":
		prompt = " [y/N] "
	else:
		raise ValueError("invalid default answer: '%s'" % default)

	while True:
		sys.stdout.write(question + prompt)
		choice = raw_input().lower()
		if default is not None and choice == '':
			return valid[default]
		elif choice in valid:
			return valid[choice]
		else:
			sys.stdout.write("Please respond with 'yes' or 'no' "\
							 "(or 'y' or 'n').\n")

def which(program):
	import os
	def is_exe(fpath):
		return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

	fpath, fname = os.path.split(program)
	if fpath:
		if is_exe(program):
			return program
	else:
		for path in os.environ["PATH"].split(os.pathsep):
			path = path.strip('"')
			exe_file = os.path.join(path, program)
			if is_exe(exe_file):
				return exe_file

	return None

def check_installation(filename,tool_name):
	if os.path.isfile(filename):
		print('%s was succesfully installed ' % tool_name)
		return True
	else:
		print 'Sorry I cannot install %s for you, install manually and try again.' % tool_name
		return False



def check_flash():
	if which('flash'):
		print '\nflash is already installed!'
		return True
	else:
		print '\nCRISPResso requires a recent version flash from: http://ccb.jhu.edu/software/FLASH/'
		if query_yes_no('Should I install flash for you?'):
			print('Ok be patient!')
			os.chdir('dependencies/')
			sb.call('tar xvzf FLASH-1.2.11.tar.gz',shell=True)
			os.chdir('FLASH-1.2.11')
			sb.call('make',shell=True)
			shutil.copy('flash', BIN_FOLDER)
			os.chdir('..')
			sb.call('rm -Rf FLASH-1.2.11',shell=True)
			os.chdir('..')
			print('flash should be installed (please check the output)')
		
	if not check_installation(os.path.join(BIN_FOLDER,'flash'),'flash'):
		sys.exit(1)
				
	else:
		return False

		
def check_needle():
	if which('needle'):
		print '\nneedle is already installed!'
		return  True
	else:
		print '\nCRISPResso requires a recent version needle from the EMBOSS suite(>=6): http://emboss.sourceforge.net/download/#Stable/'
		if query_yes_no('Should I install needle for you?'):
			print('Ok be patient!')
			os.chdir('dependencies/')
			print '\nDownloading EMBOSS source please be patient...'
			urllib.urlretrieve ("ftp://emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0.tar.gz",'EMBOSS-6.6.0.tar.gz')
			sb.call('tar xvzf EMBOSS-6.6.0.tar.gz',shell=True)
			os.chdir('EMBOSS-6.6.0')
			cmd_cfg='./configure --prefix=%s --without-x' % INSTALLATION_PATH
			print cmd_cfg    
			sb.call(cmd_cfg,shell=True)
			sb.call('make && make install',shell=True)
			os.chdir('..')
			sb.call('rm -Rf EMBOSS-6.6.0',shell=True)
			sb.call('rm EMBOSS-6.6.0.tar.gz',shell=True)
			os.chdir('..')   
			#installa needle
			print('needle should be installed (please check the output)')
		
	if not check_installation(os.path.join(BIN_FOLDER,'needle'),'needle'):
		sys.exit(1)
	else:
		return False

def install_dependencies():
		
	CURRENT_PLATFORM=platform.system().split('_')[0]

	if CURRENT_PLATFORM not in  ['Linux','Darwin'] and platform.architecture()!='64bit':
		print 'Sorry your platform is not supported\n CRISPResso is supported only on 64bit versions of Linux or OSX '
		sys.exit(1)
	
	if query_yes_no('I will install CRISPResso dependencies in:%s \n\nIs that ok?' % INSTALLATION_PATH):    
   
		if not os.path.exists(INSTALLATION_PATH):
			print 'OK, creating the folder:%s' % INSTALLATION_PATH
			os.makedirs(INSTALLATION_PATH)
			os.makedirs(BIN_FOLDER)
		else:
			print '\nI cannot create the folder!\nThe folder %s is not empty!' % INSTALLATION_PATH
			if query_yes_no('\nCan I overwrite its content? \nWARNING: all the files inside will be overwritten!'):
				#shutil.rmtree(INSTALLATION_PATH)
				#os.makedirs(INSTALLATION_PATH)
				try:
					os.makedirs(BIN_FOLDER)
				except:
					pass
			else:
				print '\nOK, install CRISPResso dependencies in a different PATH running again this script with: \n\npython setup.py YOUR_PATH'
				sys.exit(1)
		
	else:
		print '\nOK, to install CRISPResso dependencies in a different PATH just run this script again with: \n\npython setup.py YOUR_PATH'
		sys.exit(1)

	print 'CHECKING DEPENDENCIES...'
	
	flash_already_installed=check_flash()
	needle_already_installed=check_needle()
	dependencies_already_installed = ( flash_already_installed and needle_already_installed)
	
	print dependencies_already_installed
	
	if dependencies_already_installed:
	
		pass
	
	else:

		#ADD CRISPResso  dependencies to PATH
		home = expanduser("~")
		shell=os.environ["SHELL"].split('/')[-1]

		shell_profile=None
		line_to_add=None    

		if shell=='bash':
			if CURRENT_PLATFORM=='Darwin': 
				shell_profile='.bash_profile'
			else:
				shell_profile='.bashrc'
		
		elif shell=='sh' or shell=='ksh':
			shell_profile='.profile'
	
		elif shell=='tcsh':
			shell_profile='.tcshrc'

		elif shell=='csh':
			shell_profile='.cshrc'

		if shell in ['bash', 'sh','ksh']:
			line_to_add='export PATH=%s:$PATH' % BIN_FOLDER
		elif shell in ['tcsh','csh']:
			line_to_add= 'set path = ( %s $path)' % BIN_FOLDER
	
		cmd_add_path="echo '%s'  >> ~/%s" % (line_to_add,shell_profile)

		if not os.path.exists(os.path.join(home,shell_profile)) or not line_to_add in  open(os.path.join(home,shell_profile)).read():
			if query_yes_no('You need to add  %s to your PATH variable.\n\nShould I do for you?' % BIN_FOLDER):
				if shell_profile is None:
					print 'I cannot determine automatically the shell you are using. Please add the folder %s to your PATH manually!' % BIN_FOLDER
				print '\nExecuting:%s' % cmd_add_path
				sb.call(cmd_add_path,shell=True)
				sb.call(line_to_add,shell=True)
				print '\n\nINSTALLATION COMPLETED, open a NEW terminal and enjoy CRISPResso!'
			else:
				print '\nOK.\n\nNOTE: to run CRISPResso all the files in %s should be in your PATH' % BIN_FOLDER
		else:
			print '\n\nINSTALLATION COMPLETED, open a NEW terminal and enjoy CRISPResso!'    

if __name__ == '__main__':
    main()
    if sys.argv[1]=='install':
    	print 'Python package installed'
    	print '\n\nChecking dependencies...'
    	install_dependencies()
    	print '\nAll done!'









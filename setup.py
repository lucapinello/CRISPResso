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

#TO INSTALL CRISPRESSO DEPENDECIENS IN A CUSTOM LOCATION SET THE ENV VARIABLE: CRISPRESSO_DEPENDENCIES_FOLDER
if os.environ.get('CRISPRESSO_DEPENDENCIES_FOLDER'):
	INSTALLATION_PATH=os.environ.get('CRISPRESSO_DEPENDENCIES_FOLDER')
else:
	INSTALLATION_PATH='%s/CRISPResso_dependencies' % os.environ['HOME']

BIN_FOLDER=os.path.join(INSTALLATION_PATH,'bin')

def main():
	if float(sys.version[:3])<2.6 or float(sys.version[:3])>=2.8:
		sys.stdout.write("CRITICAL: Python version must be 2.7!\n")
		sys.exit(1)


	version = re.search(
    	'^__version__\s*=\s*"(.*)"',
    	open('CRISPResso/CRISPRessoCORE.py').read(),
    	re.M
    	).group(1)


	if float(sys.version[:3])<2.6 or float(sys.version[:3])>=2.8:
    		sys.stdout.write("ERROR: Python version must be 2.6 or 2.7!\n")
    		sys.exit(1)

	setup(
		  version=version,
          name = "CRISPResso",
          include_package_data = True,
    	  packages = ["CRISPResso"],
    	  package_dir={'CRISPResso': 'CRISPResso'},
          package_data={'CRISPResso': ['data/*']},
    	  entry_points = {
        	"console_scripts": ['CRISPResso = CRISPResso.CRISPRessoCORE:main',
          'CRISPRessoPooled = CRISPResso.CRISPRessoPooledCORE:main',
          'CRISPRessoWGS = CRISPResso.CRISPRessoWGSCORE:main',
          'CRISPRessoCompare = CRISPResso.CRISPRessoCompareCORE:main',
          'CRISPRessoPooledWGSCompare = CRISPResso.CRISPRessoPooledWGSCompareCORE:main',
          'CRISPRessoCount = CRISPResso.CRISPRessoCountCORE:main']
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
              'numpy>=1.9',
              'pandas>=0.15',
              'matplotlib>=1.3.1',
              'biopython>=1.6.5',
              'argparse>=1.3',
			  'seaborn>=0.7.1',
              ],

          )

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
		sys.stdout.write('%s was succesfully installed ' % tool_name)
		return True
	else:
		sys.stdout.write('Sorry I cannot install %s for you, install manually and try again.' % tool_name)
		return False



def check_flash():
	if which('flash'):
		sys.stdout.write('\nflash is already installed!')
		return True
	else:
		sys.stdout.write('\nCRISPResso requires a recent version FLASh from: http://ccb.jhu.edu/software/FLASH/')
		sys.stdout.write('\nTrying to install it, please be patient!')
		os.chdir('dependencies/')
		sb.call('tar xvzf FLASH-1.2.11.tar.gz',shell=True)
		os.chdir('FLASH-1.2.11')
		sb.call('make',shell=True)
		shutil.copy('flash', BIN_FOLDER)
		os.chdir('..')
		sb.call('rm -Rf FLASH-1.2.11',shell=True)
		os.chdir('..')
		sys.stdout.write('\nFLASh should be installed (please check the output)')

	if not check_installation(os.path.join(BIN_FOLDER,'flash'),'flash'):
		sys.exit(1)

	else:
		return False


def check_needle():
	if which('needle'):
		sys.stdout.write ('\nneedle is already installed!')
		return  True
	else:
		sys.stdout.write ('\nCRISPResso requires a recent version needle from the EMBOSS suite(>=6): http://emboss.sourceforge.net/download/#Stable/')
		sys.stdout.write('\nTrying to download and install it, please be patient!')
		os.chdir('dependencies/')
		sys.stdout.write( '\nDownloading EMBOSS source please be patient...')
		urllib.urlretrieve ("ftp://emboss.open-bio.org/pub/EMBOSS/old/6.5.0/EMBOSS-6.5.7.tar.gz",'EMBOSS-6.5.7.tar.gz')
		sb.call('tar xvzf EMBOSS-6.5.7.tar.gz',shell=True)
		os.chdir('EMBOSS-6.5.7')
		cmd_cfg='./configure --prefix=%s --without-x' % INSTALLATION_PATH
		sb.call(cmd_cfg,shell=True)
		sb.call('make && make install',shell=True)
		os.chdir('..')
		sb.call('rm -Rf EMBOSS-6.5.7',shell=True)
		sb.call('rm EMBOSS-6.5.7.tar.gz',shell=True)
		os.chdir('..')
		#installa needle
		sys.stdout.write('\nneedle should be installed (please check the output)')

	if not check_installation(os.path.join(BIN_FOLDER,'needle'),'needle'):
		sys.exit(1)
	else:
		return False

def install_dependencies():

	CURRENT_PLATFORM=platform.system().split('_')[0]

	if CURRENT_PLATFORM not in  ['Linux','Darwin'] and platform.architecture()!='64bit':
		sys.stdout.write('Sorry your platform is not supported\n CRISPResso is supported only on 64bit versions of Linux or OSX ')
		sys.exit(1)

	if not os.path.exists(INSTALLATION_PATH):
		sys.stdout.write ('OK, creating the folder:%s' % INSTALLATION_PATH)
		os.makedirs(INSTALLATION_PATH)
		os.makedirs(BIN_FOLDER)
	else:
		sys.stdout.write ('\nI cannot create the folder!\nThe folder %s is not empty!' % INSTALLATION_PATH)

		try:
			os.makedirs(BIN_FOLDER)
		except:
			pass


	sys.stdout.write( '\nCHECKING DEPENDENCIES...')

	flash_already_installed=check_flash()
	needle_already_installed=check_needle()
	dependencies_already_installed = ( flash_already_installed and needle_already_installed)


	if dependencies_already_installed:

		pass

	else:

		#ADD CRISPResso  dependencies to PATH
		home = expanduser("~")
		shell=os.environ["SHELL"].split('/')[-1]

		shell_profile=[None,]
		line_to_add=None

		if shell=='bash':
			shell_profiles=['.bash_profile','.bashrc']

		elif shell=='sh' or shell=='ksh':
			shell_profiles=['.profile']

		elif shell=='tcsh':
			shell_profiles=['.tcshrc']

		elif shell=='csh':
			shell_profiles=['.cshrc']

		if shell in ['bash', 'sh','ksh']:
			line_to_add='export PATH=%s:$PATH' % BIN_FOLDER
		elif shell in ['tcsh','csh']:
			line_to_add= 'set path = ( %s $path)' % BIN_FOLDER



		for shell_profile in shell_profiles:
			cmd_add_path="echo '%s'  >> ~/%s" % (line_to_add,shell_profile)
			if not os.path.exists(os.path.join(home,shell_profile)) or not line_to_add in  open(os.path.join(home,shell_profile)).read():

					if shell_profile is None:
						sys.stdout.write ('I cannot determine automatically the shell you are using. Please add the folder %s to your PATH manually!' % BIN_FOLDER)
					sys.stdout.write ( '\nExecuting:%s' % cmd_add_path)
					sb.call(cmd_add_path,shell=True)
					sb.call(line_to_add,shell=True)


		sys.stdout.write ('\n\nINSTALLATION COMPLETED, open a NEW terminal and enjoy CRISPResso!'    )

if __name__ == '__main__':
    main()
    if sys.argv[1]=='install':
    	sys.stdout.write ('\nPython package installed')
    	sys.stdout.write ('\n\nChecking dependencies...')
    	install_dependencies()
    	sys.stdout.write ('\nAll done!')

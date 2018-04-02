#!/usr/bin/env bash

set -e # Exit as soon as any line in the bash script fails

ROOTDIR=$( cd $(dirname $0) ; pwd -P ) # path to main PhaME directory

echo
exec &> >(tee -a  install.log)
exec 2>&1 # copies stderr onto stdout

# create a directory where all dependencies will be installed
cd $ROOTDIR
mkdir -p thirdParty
cd thirdParty
mkdir -p $ROOTDIR/bin
export PATH="$ROOTDIR/bin/":"$ROOTDIR/thirdParty/miniconda/bin/":$PATH
# Minimum Required versions of dependencies
# nucmer 3.1 is packaged in mummer 3.23
bowtie2_VER=2.2.8
bwa_VER=0.7.15
cmake_VER=3.0.1
cpanm_VER=1.7039
FastTree_VER=2.1.9
mafft_VER=7.305
miniconda_VER=4.2.12
mummer_VER=3.23
muscle_VER=3.8.31
paml_VER=4.9
R_VER=3.3.1
RAxML_VER=8.2.10
samtools_VER=1.3.1
perl5_VER=5.8.0
bcftools_VER=1.6
hyphy_VER=2.3.11
bbmap_VER=37.90

#minimum required version of Perl modules
perl_File_Basename_VER=2.85
perl_File_Path_VER=2.09
perl_Getopt_Long_VER=2.45
perl_IO_Handle_VER=1.28
perl_Parllel_ForkManager_VER=1.17
perl_Time_BaseName=2.85
perl_Time_HiRes=1.9726
perl_Statistics_Distributions_VER=1.02
perl_Test_Exception_VER=0.43



################################################################################
done_message () {
   if [ $? == 0 ]; then
      if [ "x$1" != "x" ]; then
          echo $1;
      else
          echo "done.";
      fi
   else
      echo "Installation failed." $2
      exit 1;
   fi
}

################################################################################
install_miniconda()
{
echo "--------------------------------------------------------------------------
                           downloading miniconda
--------------------------------------------------------------------------------
"

if [[ "$OSTYPE" == "darwin"* ]]
then
{

  curl -o miniconda.sh https://repo.continuum.io/miniconda/Miniconda2-4.2.12-MacOSX-x86_64.sh
  chmod +x miniconda.sh
  ./miniconda.sh -b -p $ROOTDIR/thirdParty/miniconda -f
  ln -sf $ROOTDIR/thirdParty/miniconda/bin/conda $ROOTDIR/bin/conda
  ln -sf $ROOTDIR/thirdParty/miniconda/bin/pip $ROOTDIR/bin/pip
  export PATH=$ROOTDIR/bin:$PATH
}
else
{  

  wget https://repo.continuum.io/miniconda/Miniconda2-4.2.12-Linux-x86_64.sh -O miniconda.sh
  chmod +x miniconda.sh
  ./miniconda.sh -b -p $ROOTDIR/thirdParty/miniconda -f
  ln -sf $ROOTDIR/thirdParty/miniconda/bin/conda $ROOTDIR/bin/conda
  ln -sf $ROOTDIR/thirdParty/miniconda/bin/pip $ROOTDIR/bin/pip
  export PATH=$ROOTDIR/bin:$PATH

}
fi
}

################################################################################
checkSystemInstallation()
{
    IFS=:
    for d in $PATH; do
      if test -x "$d/$1"; then return 0; fi
    done
    return 1
}


checkLocalInstallation()
{
    IFS=:
    for d in $ROOTDIR/thirdParty/miniconda/bin; do
      if test -x "$d/$1"; then return 0; fi
    done
    return 1
}

checkHyPhyInstallation()
{
    IFS=:
    for d in $ROOTDIR/bin; do
      if test -x "$d/$1"; then return 0; fi
    done
    return 1
}

checkPerlModule()
{
   perl -e "use lib \"$ROOTDIR/lib\"; use $1;"
   return $?
}

################################################################################
install_cpan()
{
echo "--------------------------------------------------------------------------
                           installing cpanm v2.00
--------------------------------------------------------------------------------
"
conda install --yes c dan_blanchard perl-cpan=2.00 -p $ROOTDIR/thirdParty/miniconda
ln -sf $ROOTDIR/thirdParty/miniconda/bin/cpan $ROOTDIR/bin/cpan
echo "--------------------------------------------------------------------------
                           installed cpanm v2.00
--------------------------------------------------------------------------------
"

}

################################################################################
install_cpanm()
{
echo "--------------------------------------------------------------------------
                           installing cpanm v1.7039 
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda perl-app-cpanminus=$cpanm_VER -p $ROOTDIR/thirdParty/miniconda
ln -sf $ROOTDIR/thirdParty/miniconda/bin/cpanm $ROOTDIR/bin/cpanm
echo "--------------------------------------------------------------------------
                          cpanm v1.7039 installed
--------------------------------------------------------------------------------
"
}

install_mummer()
{
echo "--------------------------------------------------------------------------
                           installing mummer v3.23 or nucmer 3.21
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda mummer=$mummer_VER -p $ROOTDIR/thirdParty/miniconda
ln -sf $ROOTDIR/thirdParty/miniconda/opt/mummer-$mummer_VER/nucmer $ROOTDIR/bin/nucmer
ln -sf $ROOTDIR/thirdParty/miniconda/bin/mummer $ROOTDIR/bin/mummer
ln -sf $ROOTDIR/thirdParty/miniconda/bin/delta-filter $ROOTDIR/bin/delta-filter
ln -sf $ROOTDIR/thirdParty/miniconda/bin/show-snps $ROOTDIR/bin/show-snps
ln -sf $ROOTDIR/thirdParty/miniconda/bin/show-coords $ROOTDIR/bin/show-coords
# conda install --yes -c bioconda mummer=$mummer_VER -p $ROOTDIR/thirdParty/miniconda
echo "
------------------------------------------------------------------------------
                           mummer v3.23 installed or nucmer 3.21
------------------------------------------------------------------------------
"
}
install_bbmap()
{
echo "--------------------------------------------------------------------------
                           installing bbmap v37.90
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda bbmap -p $ROOTDIR/thirdParty/miniconda
ln -sf $ROOTDIR/thirdParty/miniconda/bin/bbmap.sh $ROOTDIR/bin/bbmap.sh
ln -sf $ROOTDIR/thirdParty/miniconda/bin/sketch.sh $ROOTDIR/bin/sketch.sh
ln -sf $ROOTDIR/thirdParty/miniconda/bin/comparesketch.sh $ROOTDIR/bin/comparesketch.sh

echo "
------------------------------------------------------------------------------
                           mummer v37.90 
------------------------------------------------------------------------------
"
}


install_cmake()
{
echo "--------------------------------------------------------------------------
                           installing cmake v3.0.1
--------------------------------------------------------------------------------
"
conda install --yes -c anaconda cmake=$cmake_VER -p $ROOTDIR/thirdParty/miniconda
ln -sf $ROOTDIR/thirdParty/miniconda/bin/cmake $ROOTDIR/bin/cmake
echo "
------------------------------------------------------------------------------
                           cmake v3.0.1 installed
------------------------------------------------------------------------------
"
}

install_bwa()
{
echo "------------------------------------------------------------------------------
                           Downloading bwa v0.7.15
------------------------------------------------------------------------------
"
conda install --yes -c bioconda bwa=$bwa_VER -p $ROOTDIR/thirdParty/miniconda
ln -sf $ROOTDIR/thirdParty/miniconda/bin/bwa $ROOTDIR/bin/bwa
echo "
------------------------------------------------------------------------------
                           bwa v0.7.15 installed
------------------------------------------------------------------------------
"
}

install_bowtie2()
{
echo "--------------------------------------------------------------------------
                           installing bowtie2 v2.2.8
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda bowtie2=$bowtie2_VER -p $ROOTDIR/thirdParty/miniconda
ln -sf $ROOTDIR/thirdParty/miniconda/bin/bowtie2 $ROOTDIR/bin/bowtie2
echo "
------------------------------------------------------------------------------
                           bowtie2 v2.2.8 installed
------------------------------------------------------------------------------
"
}

install_samtools()
{
echo "--------------------------------------------------------------------------
                           Downloading samtools v1.3.1
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda samtools=$samtools_VER -p $ROOTDIR/thirdParty/miniconda
ln -sf $ROOTDIR/thirdParty/miniconda/bin/samtools $ROOTDIR/bin/samtools
echo "
--------------------------------------------------------------------------------
                           samtools v1.3.1 installed
--------------------------------------------------------------------------------
"
}

install_FastTree()
{
echo "--------------------------------------------------------------------------
                           Downloading FastTree v2.1.9
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda FastTree=$FastTree_VER -p $ROOTDIR/thirdParty/miniconda
ln -sf $ROOTDIR/thirdParty/miniconda/bin/FastTree $ROOTDIR/bin/FastTree
ln -sf $ROOTDIR/thirdParty/miniconda/bin/FastTreeMP $ROOTDIR/bin/FastTreeMP
echo "
--------------------------------------------------------------------------------
                           FastTree v2.1.9
--------------------------------------------------------------------------------
"
}

install_RAxML()
{
echo "--------------------------------------------------------------------------
                           Downloading RAxML v$RAxML_VER
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda RAxML=$RAxML_VER -p $ROOTDIR/thirdParty/miniconda
ln -sf $ROOTDIR/thirdParty/miniconda/bin/raxmlHPC $ROOTDIR/bin/raxmlHPC
ln -sf $ROOTDIR/thirdParty/miniconda/bin/raxmlHPC $ROOTDIR/bin/raxmlHPC-PTHREADS
echo "
--------------------------------------------------------------------------------
                           RAxML v$RAxML_VER
--------------------------------------------------------------------------------
"
}

install_muscle()
{
echo "--------------------------------------------------------------------------
                           Downloading muscle v3.8.31
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda muscle=$muscle_VER -p $ROOTDIR/thirdParty/miniconda
ln -sf $ROOTDIR/thirdParty/miniconda/bin/muscle $ROOTDIR/bin/muscle
echo "
--------------------------------------------------------------------------------
                           muscle v3.8.31
--------------------------------------------------------------------------------
"
}

install_mafft()
{
echo "--------------------------------------------------------------------------
                           Downloading mafft v7.305
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda mafft=$mafft_VER -p $ROOTDIR/thirdParty/miniconda
ln -sf $ROOTDIR/thirdParty/miniconda/bin/mafft $ROOTDIR/bin/mafft
echo "
--------------------------------------------------------------------------------
                           mafft v7.305
--------------------------------------------------------------------------------
"
}

install_paml()
{
echo "--------------------------------------------------------------------------
                           Downloading paml v4.9
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda paml=$paml_VER -p $ROOTDIR/thirdParty/miniconda
echo "
--------------------------------------------------------------------------------
                           paml v4.9
--------------------------------------------------------------------------------
"
}

install_bcftools()
{
echo "--------------------------------------------------------------------------
                           installing bcftools v$bcftools_VER
--------------------------------------------------------------------------------
"
conda install --yes -c bioconda bcftools=$bcftools_VER -p $ROOTDIR/thirdParty/miniconda
ln -sf $ROOTDIR/thirdParty/miniconda/bin/bcftools $ROOTDIR/bin/bcftools
echo "--------------------------------------------------------------------------
                           installed bcftools v$bcftools_VER
--------------------------------------------------------------------------------
"
}

install_hyphy()
{
echo "--------------------------------------------------------------------------
                           installing hyphy v$hyphy_VER
--------------------------------------------------------------------------------
"
#curl -o thirdParty/hyphy-2.3.11.zip https://github.com/veg/hyphy/archive/2.3.11.zip
cd $ROOTDIR
unzip thirdParty/hyphy-2.3.11.zip -d thirdParty/;
gcc -v;
which gcc;
cd $ROOTDIR/thirdParty/hyphy-2.3.11;
cmake -DINSTALL_PREFIX=$ROOTDIR/thirdParty/hyphy-2.3.11
make MP
make install
make GTEST
cd $ROOTDIR
ln -sf $ROOTDIR/thirdParty/hyphy-2.3.11/HYPHYMP $ROOTDIR/bin/HYPHYMP
echo "--------------------------------------------------------------------------
                           installed hyphy v$hyphy_VER
--------------------------------------------------------------------------------
"
}

install_perl_Getopt_Long()
{
echo "--------------------------------------------------------------------------
                installing Perl Module Getopt::Long v2.45
--------------------------------------------------------------------------------
"
cpanm Getopt::Long@$perl_Getopt_Long_VER
echo "
--------------------------------------------------------------------------------
                           Getopt::Long v2.45 installed
--------------------------------------------------------------------------------
"
}
install_perl_Test_Exception()
{
echo "--------------------------------------------------------------------------
                installing Perl Module Test::Exception v0.46
--------------------------------------------------------------------------------
"
cpanm Test::Exception@$perl_Test_Exception_VER
echo "
--------------------------------------------------------------------------------
                           Test::Exception v0.46 installed
--------------------------------------------------------------------------------
"
}

install_perl_Statistics_Distributions()
{
echo "--------------------------------------------------------------------------
          installing Perl Module Statistics::Distributions v1.02
--------------------------------------------------------------------------------
"
cpanm Statistics::Distributions@$perl_Statistics_Distributions_VER
echo "
--------------------------------------------------------------------------------
                          Statistics::Distributions v1.02 installed
--------------------------------------------------------------------------------
"
}

install_perl_Time_HiRes()
{
echo "--------------------------------------------------------------------------
                           installing Perl Module Time::HiRes v1.9726
--------------------------------------------------------------------------------
"
cpanm Time::HiRes@$perl_Time_HiRes_VER
echo "
--------------------------------------------------------------------------------
                          Time::HiRes v1.9726 installed
--------------------------------------------------------------------------------
"
}

install_perl_File_Path()
{
echo "--------------------------------------------------------------------------
                           installing Perl Module File::Path v2.09
--------------------------------------------------------------------------------
"
cpanm File::Path@$perl_File_Path_VER
echo "
--------------------------------------------------------------------------------
                          File::Path v2.09 installed 
--------------------------------------------------------------------------------
"
}

install_perl_File_Basename()
{
echo "--------------------------------------------------------------------------
                           installing Perl Module File::Basename v2.85
--------------------------------------------------------------------------------
"
cpanm File::Basename@$perl_File_Basename_VER
echo "
--------------------------------------------------------------------------------
                          File::Basename Resinstalled v2.85
--------------------------------------------------------------------------------
"
}

install_perl_File_Copy()
{
echo "--------------------------------------------------------------------------
                           installing Perl Module File::Copy v2.30
--------------------------------------------------------------------------------
"
cpanm File::Copy@$perl_File_Copy_VER
echo "
--------------------------------------------------------------------------------
                          File::Copy v2.30 installed 
--------------------------------------------------------------------------------
"
}

install_perl_IO_Handle()
{
echo "--------------------------------------------------------------------------
                       installing Perl Module IO::Handle v1.35
--------------------------------------------------------------------------------
"
cpanm IO::Handle@$perl_IO_Handle_VER
echo "
--------------------------------------------------------------------------------
                           IO::Handle v1.35 installed
--------------------------------------------------------------------------------
"
}


install_perl_Parallel_ForkManager()
{
echo "--------------------------------------------------------------------------
                           installing Perl Module Parallel::ForkManager v1.17
--------------------------------------------------------------------------------
"
cpanm Parallel::ForkManager@$perl_Parllel_ForkManager_VER -l $ROOTDIR/ext
echo "
--------------------------------------------------------------------------------
                           Parallel::ForkManager v1.17 installed
--------------------------------------------------------------------------------
"
}
################################################################################
if ( checkSystemInstallation conda )
then
  echo "conda is found"
  if [ -d "$ROOTDIR/thirdParty/miniconda" ]; then
    echo "conda is installed and pointed to right environment"  
  else
    echo "Creating a separate conda enviroment ..."
    conda create --yes -p $ROOTDIR/thirdParty/miniconda
  fi
else
  echo "conda was not found"
  install_miniconda
fi

################################################################################
if ( checkSystemInstallation perl )
then
  perl_installed_VER=`perl -v 2>&1| grep 'version' | perl -nle 'print $& if m{\d+\.\d+\.\d+}'`;
  if  ( echo $perl_installed_VER $perl_VER | awk '{if($1>=$2) exit 0; else exit 1}' )
  then 
    echo " - found perl $perl_installed_VER"
  else
  echo "Required version of perl $perl_VER was not found"
  install_perl
  fi
else 
  echo "perl was not found"
  install_perl
fi

################################################################################
if ( checkSystemInstallation cpan )
then
    echo " - found cpan"
else
  echo "cpan was found"
  install_cpan
fi
################################################################################
if ( checkSystemInstallation cpanm )
then
  cpanm_installed_VER=`cpanm -V 2>&1| head -n 1 | grep 'version' | perl -nle 'print $& if m{version \d+\.\d+}'`;
  if  ( echo $cpanm_installed_VER $cpanm_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
  then 
    echo " - found cpanm $cpanm_installed_VER"
  else
  echo "Required version of cpanm was not found"
  install_cpanm
  fi
else 
  echo "cpanm was not found"
  install_cpanm
fi

################################################################################
if ( checkSystemInstallation cmake )
then
  cmake_installed_VER=`cmake -V 2>&1| head -n 1 | grep 'version' | perl -nle 'print $& if m{version \d+\.\d+}'`;
  if  ( echo $cmake_installed_VER $cmake_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
  then 
    echo " - found cmake $cmake_installed_VER"
  else
  echo "Required version of cmake was not found"
  install_cmake
  fi
else 
  echo "cmake was not found"
  install_cmake
fi

################################################################################
if ( checkSystemInstallation bbmap.sh )
then
  bbmap_installed_VER=`bbmap.sh -v 2>&1| grep 'version' | perl -nle 'print $& if m{version \d+\.\d+}'`;
  if  ( echo $bbmap_installed_VER $bbmap_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
  then 
    echo " - found bbmap $bbmap_installed_VER"
  else
  echo "Required version of bbmap was not found"
  install_bbmap
  fi
else 
  echo "bbmap was not found"
  install_bbmap
fi
################################################################################
if ( checkSystemInstallation nucmer )
then
  mummer_installed_VER=`nucmer -v 2>&1| grep 'version' | perl -nle 'print $& if m{version \d+\.\d+}'`;
  if  ( echo $mummer_installed_VER $mummer_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
  then 
    echo " - found nucmer $mummer_installed_VER"
  else
  echo "Required version of nucmer was not found"
  install_mummer
  fi
else 
  echo "nucmer was not found"
  install_mummer
fi
################################################################################

if ( checkSystemInstallation bwa )
then
bwa_installed_VER=`bwa 2>&1| grep 'Version'  | perl -nle 'print $& if m{Version: \d+\.\d+\.\d+}'`;
  if  ( echo $bwa_installed_VER $bwa_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
  then
    echo " - found BWA $bwa_installed_VER"
  else
    install_bwa
  fi
else
  echo "bwa is not found"
  install_bwa
fi
################################################################################

if ( checkSystemInstallation bowtie2 )
then
bowtie2_installed_VER=`bowtie2 --version 2>&1 | grep -v 'Compiler' |grep 'version' | perl -nle 'print $& if m{version \d+\.\d+\.\d+}'`;
  if (echo $bowtie2_installed_VER $bowtie2_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
  then
    echo " - found bowtie2 $bowtie2_installed_VER"
  else
    echo "Required version of bowtie2 $bowtie2_VER was not found"
    install_bowtie2
  fi
else
  echo "bowtie2 is not found"
  install_bowtie2
fi
################################################################################
if ( checkSystemInstallation samtools )
then
  samtools_installed_VER=`samtools 2>&1| grep 'Version'|perl -nle 'print $& if m{Version: \d+\.\d+.\d+}'`;
    if ( echo $samtools_installed_VER $samtools_VER| awk '{if($2>=$3) exit 0; else exit 1}' )
    then
    echo " - found samtools $samtools_installed_VER"
    else
    echo "Required version of samtools $samtools_VER was not found"
    install_samtools
    fi
else
  install_samtools
fi
################################################################################
if ( checkSystemInstallation FastTree )
then
  FastTree_installed_VER=`FastTree 2>&1| grep 'version'|perl -nle 'print $& if m{version \d+\.\d+.\d+}'`;
  if ( echo $FastTree_installed_VER $FastTree_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
  then
    echo " - found FastTree $FastTree_VER"
  else
    echo "Required version of $FastTree_VER was not found"
    install_FastTree
  fi
  FastTreeMP_installed_VER=`FastTreeMP 2>&1| grep 'version'|perl -nle 'print $& if m{version \d+\.\d+.\d+}'`;
  if ( echo $FastTreeMP_installed_VER $FastTree_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
  then
    echo " - found FastTreeMP $FastTreeMP_VER"
  else
    echo "Required version of FastTreeMP $FastTree_VER was not found"
    install_FastTree
  fi
else
    echo "FastTree was not found"
    install_FastTree
fi

################################################################################
if ( checkSystemInstallation raxmlHPC )
then
  RAxMLHPC_installed_VER=`raxmlHPC-PTHREADS -version 2>&1 | grep 'version' | perl -nle 'print $& if m{version \d+\.\d+\.\d+}'`
  if ( echo $RAxMLHPC_installed_VER $RAxML_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
  then 
    echo " - found RAxML $RAxMLHPC_installed_VER"
  else
    echo "Required version of $RAxML_VER was not found"
      install_RAxML
  fi  
else
  echo "RAxML was not found"
  install_RAxML
fi

################################################################################
if ( checkSystemInstallation muscle )
then
muscle_installed_VER=`muscle -version 2>&1 | grep 'v' | perl -nle 'print $& if m{v\d+\.\d+\.\d+}'`;
  if ( echo $muscle_installed_VER $muscle_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
  then
    echo " - found muscle $muscle_installed_VER"
  else
    echo "Required version of $muscle_VER was not found"
    install_muscle
  fi
else
  echo "muscle was not found"
  install_muscle
fi

################################################################################
if ( checkSystemInstallation mafft )
then
  mafft_installed_VER=`mafft --version 2>&1 | grep 'v' | perl -nle 'print $& if m{v\d+\.\d+}'`;
  if ( echo $mafft_installed_VER $mafft_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
  then 
    echo " - found mafft $mafft_installed_VER"  
  else
  echo "Required version of mafft $mafft_VER was not found"
  install_mafft
  fi
else
  echo "mafft was not found"
  install_mafft
fi


################################################################################
if ( checkSystemInstallation bcftools )
then
  bcftools_installed_VER=`bcftools -v 2>&1| grep 'bcftools' | perl -nle 'print $& if m{\d+\.\d+}'`;
  if  ( echo $bcftools_installed_VER $bcftools_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
  then 
    echo " - found bcftools $bcftools_installed_VER"
  else
  echo "Required version of bcftools was not found"
  install_bcftools
  fi
else 
  echo "bcftools was not found"
  install_bcftools
fi
################################################################################
if ( checkSystemInstallation codeml )
then
  paml_installed_VER=`evolver \0 2>&1 | grep 'version' | perl -nle 'print $& if m{version \d+\.\d+}'`;
  if ( echo $paml_installed_VER $paml_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
  then 
    echo " - found paml $paml_installed_VER"
  else 
    echo "Required version of paml $paml_VER was not found"
    install_paml
  fi
else
  echo "paml is not found"
  install_paml
fi

################################################################################
if ( checkHyPhyInstallation HYPHYMP )
then
  hyphy_installed_VER=`thirdParty/hyphy-2.3.11/HYPHYMP -v 2>&1 | grep 'HYPHY' | perl -nle 'print $& if m{HYPHY \d+\.\d+.\d\d}'`;
  if ( echo $hyphy_installed_VER $hyphy_VER | awk '{if($2>=$3) exit 0; else exit 1}' )
  then 
    echo " - found HyPhy $hyphy_installed_VER"
  else 
    echo "Required version of HyPhy $hyphy_VER was not found"
    install_hyphy
  fi
else
  echo "HyPhy was not found"
  install_hyphy
fi

################################################################################
#                                 Perl modules
################################################################################

if ( checkPerlModule Getopt::Long )
then
  echo " - found Perl module Getopt::Long"
  perl_Getopt_Long_installed_VER=`perl -MGetopt::Long -e 'print $Getopt::Long::VERSION ."\n";'`
  if ( echo $perl_Getopt_Long_installed_VER $perl_Getopt_Long_VER | awk '{if($2>=$3) exit 0; else exit 1}')
  then
    echo " - found Perl module Getopt::Long $perl_Getopt_Long_installed_VER"
  else
    echo "Required version of Getopt::Long $perl_Getopt_Long_VER was not found"
    install_perl_Getopt_Long
  fi
else
  echo "Perl Getopt::Long was not found"
  install_perl_Getopt_Long
fi

#------------------------------------------------------------------------------#

if ( checkPerlModule Time::HiRes )
then
  perl_Time_HiRes_installed_VER=`perl -MTime::HiRes -e 'print $Time::HiRes::VERSION ."\n";'`
  if ( echo $perl_Time_HiRes_installed_VER $perl_Time_HiReS_VER | awk '{if($2>=$3) exit 0; else exit 1}')
  then
    echo " - found Perl module Time::HiRes $perl_Time_HiRes_installed_VER"
  else
    echo "Required version of Time::HiRes $perl_Time_HiRes_VER was not found"
    install_perl_Time_HiRes
  fi
else
  echo "Perl Time::HiRes was not found"
  install_perl_Time_HiRes
fi
#------------------------------------------------------------------------------#
if ( checkPerlModule File::Path )
then
  perl_File_Path_installed_VER=`perl -MFile::Path -e 'print $File::Path::VERSION ."\n";'`
  if ( echo $perl_File_Path_installed_VER $perl_File_Path_VER | awk '{if($2>=$3) exit 0; else exit 1}')
  then
    echo " - found Perl module File::Path $perl_File_Path_installed_VER"
  else
    echo "Required version of File::Path $perl_File_Path_VER was not found"
    install_perl_File_Path
  fi
else
  echo "Perl File::Path is not found"
  install_perl_File_Path
fi
#------------------------------------------------------------------------------#
if ( checkPerlModule File::Basename )
then
  perl_File_BaseName_installed_VER=`perl -MFile::Basename -e 'print $File::Basename::VERSION ."\n";'`
  if ( echo $perl_File_BaseName_installed_VER $perl_File_BaseName_VER | awk '{if($2>=$3) exit 0; else exit 1}')
  then
    echo " - found Perl module File::BaseName $perl_File_BaseName_installed_VER"
  else
    echo "Required version of File::BaseName $perl_File_BaseName_VER was not found"
    install_perl_File_BaseName
  fi
else
  echo "Perl File::Basename is not found"
  install_perl_File_Basename
fi
#------------------------------------------------------------------------------#
if ( checkPerlModule File::Copy )
then
  perl_File_Copy_installed_VER=`perl -MFile::Copy -e 'print $File::Copy::VERSION ."\n";'`
  if ( echo $perl_File_Copy_installed_VER $perl_File_Copy_VER | awk '{if($2>=$3) exit 0; else exit 1}')
  then
    echo " - found Perl module File::Copy $perl_File_Copy_installed_VER"
  else
    echo "Required version of File::Copy $perl_File_Copy_VER was not found"
    install_perl_File_Copy
  fi
else
  echo "Perl File::Copy is not found"
  install_perl_File_Copy
fi
#------------------------------------------------------------------------------#
if ( checkPerlModule IO::Handle )
then
  #perl -MIO::Handle -e 'print $IO::Handle::VERSION ."\n";'
  perl_IO_Handle_installed_VER=`perl -MIO::Handle -e 'print $IO::Handle::VERSION ."\n";'`
  if ( echo $perl_IO_Handle_installed_VER $perl_IO_Handle_VER | awk '{if($2>=$3) exit 0; else exit 1}')
  then
    echo " - found Perl module IO::Handle $perl_IO_Handle_installed_VER"
  else
    echo "Required version of IO::Handle $perl_IO_Handle_VER was not found"
    install_perl_IO_Handle
  fi
else
  echo "Perl IO::Handle is not found"
  install_perl_IO_Handle
fi
#------------------------------------------------------------------------------#
if ( checkPerlModule Test::Exception )
then
  #perl -MTest::Exception -e 'print $Test::Exception::VERSION ."\n";'
  perl_Test_Exception_installed_VER=`perl -MTest::Exception -e 'print $Test::Exception::VERSION ."\n";'`
  if ( echo $perl_Test_Exception_installed_VER $perl_Test_Exception_VER | awk '{if($2>=$3) exit 0; else exit 1}')
  then
    echo " - found Perl module Test::Exception $perl_Test_Exception_installed_VER"
  else
    echo "Required version of Test::Exception $perl_Test_Exception_VER was not found"
    install_perl_Test_Exception
  fi
else
  echo "Perl Test::Exception is not found"
  install_perl_Test_Exception
fi
#------------------------------------------------------------------------------#
if ( checkPerlModule Parallel::ForkManager )
then
  perl_Parallel_ForkManager_installed_VER=`perl -MParallel::ForkManager -e 'print $Parallel::ForkManager::VERSION ."\n";'`
  if ( echo $perl_Parallel_ForkManager_installed_VER $perl_Parallel_ForkManager_VER | awk '{if($1>=$2) exit 0; else exit 1}')
  then
    echo " - found Perl module Parallel::ForkManager $perl_Parallel_ForkManager_installed_VER"
  else
    echo "Required version of Parallel::ForkManager $perl_Parallel_ForkManager_VER was not found"
    install_perl_Parallel_ForkManager
  fi
else
  echo "Perl Parallel::ForkManager is not found"
  install_perl_Parallel_ForkManager
fi
#------------------------------------------------------------------------------#
if ( checkPerlModule Statistics::Distributions )
then
  perl_Statistics_Distributions_installed_VER=`cpan -D Statistics::Distributions 2>&1 | grep 'Installed' | perl -nle 'print $& if m{Installed: \d+\.\d+}'`
  if ( echo $perl_Statistics_Distributions_installed_VER $perl_Statistics_Distributions_VER | awk '{if($2>=$3) exit 0; else exit 1}')
  then
    echo " - found Perl module Statistics::Distributions $perl_Statistics_Distributions_installed_VER"
  else
    echo "Required version of Statistics::Distributions $perl_Statistics_Distributions_VER was not found"
    install_perl_Statistics_Distributions
  fi
else
  echo "Perl Statistics::Distributions is not found"
  install_perl_Statistics_Distributions
fi

################################################################################




################################################################################
echo "
================================================================================
                 PhaME installed successfully.
================================================================================
Check phame.ctl for the control file

Quick start:
    perl src/runPhaME.pl test/phame.ctl
    Please use the full path to specify the location of fasta file.
    Read README.md for detail
Check our github site for update:
    https://github.com/LANL-Bioinformatics/PhaME
";

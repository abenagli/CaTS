use Env;
use Cwd;



open(USERCONFIG,$ARGV[0]);
while(<USERCONFIG>)
{
  chomp;
  s/#.*//;                # no comments
  s/^\s+//;               # no leading white
  s/\s+$//;               # no trailing white
  my($var,$value) = split(/\s*=\s*/, $_, 2);
  $User_Preferences{$var} = $value;
}

$LISTFile    = $User_Preferences{"LISTFile"};
$GDMLFile    = $User_Preferences{"GDMLFile"};
$MACFile     = $User_Preferences{"MACFile"};
$OUTPUTFolder= $User_Preferences{"OUTPUTFolder"};
$OUTPUTLabel = $User_Preferences{"OUTPUTLabel"};
$EXEName     = $User_Preferences{"EXEName"};

$ISLxplus   = $User_Preferences{"ISLxplus"};
$ISHercules = $User_Preferences{"ISHercules"};
$QUEUE      = $User_Preferences{"QUEUE"};



open(LANCIA,">","lancia.sh") or die "Can't open file lancia.sh";
print LANCIA "#! /bin/zsh \n";



$CURRENTDir = getcwd();

$OUTPUTDir = $OUTPUTFolder."/".$OUTPUTLabel."/";
if( $ISLxplus )
{
  print("/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select mkdir -p /eos/cms/".$OUTPUTDir."\n");
  system("/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select mkdir -p /eos/cms/".$OUTPUTDir."/");
}

$OUTPUTJobDir = "";
if( $ISLxplus )
{
  $OUTPUTJobDir = $CURRENTDir."/jobs/".$OUTPUTLabel;
  if( ! -e $OUTPUTJobDir )
  {
    system("mkdir ".$OUTPUTJobDir);
  }
}



open(LIST,$LISTFile);
while(<LIST>)
{
  chomp;
  s/#.*//;                # no comments
  s/^\s+//;               # no leading white
  s/\s+$//;               # no trailing white
  
  ($label,$particle,$energy,$Nevts,$Njobs,$Nfirst) = split(" ");
  
  if( $label eq "" )
  {
    next;
  }
  
  print($label."_".$particle," ".$energy." ".$Nevts." ".$Njobs." ".$Nfirst."\n");
  
  
  
  $STAGEOUTDir = $OUTPUTDir."/run_".$label;
  if( $ISLxplus )
  {
    print("/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select mkdir -p /eos/cms/".$STAGEOUTDir."/\n");
    system("/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select mkdir -p /eos/cms/".$STAGEOUTDir."/");
  }
  
  
  
  $runDir = $OUTPUTJobDir."/run_".$label."/";
  if( ! -e $runDir )
  {
    #print("mkdir ".$runDir."\n");
    system("mkdir ".$runDir);
  }
  
  for($iJob = $Nfirst; $iJob < $Nfirst+$Njobs; ++$iJob)
  {
    $jobDir = $runDir."/job_".$iJob."/";
    if( ! -e $jobDir )
    {
      #print("mkdir ".$jobDir."\n");
      system("mkdir ".$jobDir);
    }
    
    else
    {
      if( $ARGV[1] == 1 )
      {
        #print("rm -rf ".$jobDir."\n");
        system("rm -rf ".$jobDir);
        #print("mkdir ".$jobDir."\n");
        system("mkdir ".$jobDir);
      }
      else
      {
        print(">>> output directory not empty. Remove it or run ./launchJobs.pl params.cfg 1 <<< \n");
        die;
      }
    }
    
    
    $jobMacFile = $jobDir."/particle.mac";
    system("cat ".$MACFile." | sed -e s%PARTICLE%".$particle.
                         "%g | sed -e s%ENERGY%".$energy.
                         "%g | sed -e s%NEVTS%".$Nevts.
                         "%g | sed -e s%JOB%".$iJob.
                         "%g > ".$jobMacFile);
    $jobOut = $jobDir."/log_".$label."_".$iJob.".txt";
    
    
    $jobSh = $jobDir."/job_".$label."_".$iJob.".sh";
    open(JOBSH,">",$jobSh) or die "Can't open file ".$jobSh;
    
    print JOBSH "echo \$SHELL\n";
    print JOBSH "pwd \n";
    print JOBSH "touch ".$jobDir."/job.start\n";
    
    if( $ISLxplus )
    {
      print JOBSH "source /afs/cern.ch/sw/lcg/external/gcc/4.7/x86_64-slc6-gcc47-opt/setup.sh\n";
      print JOBSH ". /afs/cern.ch/sw/lcg/external/geant4/10.0.p03/x86_64-slc6-gcc47-opt/CMake-setup.sh\n";
      print JOBSH "touch ".$jobDir."/job.run\n";
      print JOBSH "unbuffer ".$EXEName." ".$GDMLFile." ".$jobMacFile."\n";
      print JOBSH "/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select cp ".$particle."_".$energy."GeV_hits_".$iJob.".root /eos/cms/".$STAGEOUTDir."/\n";
      print JOBSH "touch ".$jobDir."/job.done\n";
    }
    
    if( $ISHercules == 1 )
    {
      print JOBSH "cd ".$jobDir."\n";
      print JOBSH "pwd \n";
      print JOBSH "touch ".$jobDir."/job.run\n";
      print JOBSH "unbuffer ".$EXEName." ".$jobDir."job.cfg out_".$iJob." >& ".$jobOut."\n";
      print JOBSH "touch ".$jobDir."/job.done\n";
    }
    
    system("chmod +x ".$jobSh);
    
    if( $ISLxplus == 1 )
    {
      print LANCIA "bsub -cwd ".$jobDir." -q ".$QUEUE." ".$jobSh."\n";
    }
  }
}

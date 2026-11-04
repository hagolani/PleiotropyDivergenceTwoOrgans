#!/bin/sh

# SLURM options:

#SBATCH --job-name=evox    # Job name
#SBATCH --output=evox.log   # Standard output and error log

#SBATCH --partition=E5               # Partition choice (htc by default)
#SBATCH --nodes=1

#SBATCH --ntasks=16                    # Run a single task
#SBATCH --mem=10000                    # Memory in MB per default
#SBATCH --time=7-00:00:00             # Max time limit = 7 days

#SBATCH --mail-user=pfhagolani@gmail.com   # Where to send mail
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)

# Commands to be submitted:


#To modify
maxCPU=16  #how many cores will you use

psize=200  #population size
mut=0.025 #0.025  #mutation chance of each parameter (pi)
ngen=200    #number of generations
mutstep=7.5 #multiplies the step matrix for mutations. Affects all parameters equally
pleio=1   #do teeth change independently? 1 independent ; else mutations the same for both


#if hox gene is present it will still be used, that line needs to be taken out if you do not want to use it in developing the teeth (use "usehox" below)

hoxmutt=0.1 #chance modifier to mutate hox vector (hi) (chance is hoxmutt*mut)
hoxhnmut=0.05 #chance modifier to mutate hox gene (hn) (chance is hoxhnmut*mut) ; about 1 every two generations

echo $hoxmutt $hoxhnmut

fittype=2 #fitness used 1 distance or 2 surface

hox=1 #mutate hox gene? 1 yes, else no. If trying to use with input file that does not have the line, it will simply run as if no hox gene
usehox1=1 #if 1 we add hox if not already present in file. If ==2 we take it out. For teeth 1
usehox2=1 #if 1 we add hox if not already present in file. If ==2 we take it out. For teeth 2

hoxstep=1  #multiplies mutation step for hox. step is = step*mutstep*hoxstep

conserve=0

#Making file with evo details
rm -f evo_parameters.dat
touch evo_parameters.dat
echo "Population_size $psize" >> evo_parameters.dat
echo "Mutation_chance_per_parameter $mut" >> evo_parameters.dat
echo "Mutation_step_multiplier $mutstep" >> evo_parameters.dat
echo "Generations $ngen" >> evo_parameters.dat
#echo -n "Parentals " >> evo_parameters.dat
echo "Parentals $(cat initial_file_name)" >> evo_parameters.dat
echo "Pleiotropic_1yes_elseNo $pleio" >> evo_parameters.dat 
echo "mutate_hox_gene_1yes_elseNo $hox" >> evo_parameters.dat 
echo "modifier_mutate_hox_vector $hoxmutt" >> evo_parameters.dat
echo "modifier_mutate_hox_hn $hoxhnmut" >> evo_parameters.dat
echo "fitness_used_1distance_2surface $fittype" >> evo_parameters.dat
echo "use_hox_gene_1yes_elseNo_forEachTeeth $usehox1 $usehox2" >> evo_parameters.dat 
echo "conservative_evo_until_generation $conserve" >> evo_parameters.dat 

#Please do not touch
one=1
two=2
#model stuff, but should not change unless you have a reason
ini=0 #for 1 genome population introduces a bit of "noise" in the parameters, for 2 genomes if makes mutation equal to 0 just to make a copy of the 1 genome population
half=$((psize / two))  
psize_run=$((psize * two))  #independent teeth, so each "individual/mouth" develops two teeth

maxCPU2=$((maxCPU * two)) #timeout of each runt.e included in "ps -aef | grep runt.e | wc -l"


###########

fittype2=$fittype  #conservative, until generation == conserve


#make initial list
#TODO

#echo "sealmaster1.txt sealmaster2.txt" > "TOMUT.txt"
#psize00=`expr $psize - $one`
#for s in $(seq 1 $psize00); do   
#   echo "sealmaster1.txt sealmaster2.txt" >> "TOMUT.txt" #each individual in a row, ie, two teeth per row
#done

# 1. Make the initial population, should be all very similar, all parameters with very small changes in respect to the original sealmaster.txt
# Then runs the rest of generations

for gen in $(seq 0 $ngen); do #number of generations
  echo "Generation" $gen
  i=0
  psize0=$(wc -l TOMUT.txt | awk '{ print $1 }') #number of files to pontentially mutate 
  echo $psize0
  while [ $i -lt $psize0 ]
   do
       i=`expr $i + $one`
       name=$(sed -n "$i"p TOMUT.txt) #name includes both teeth
        echo $name $ini $mut $gen $i $mutstep $pleio $hox $hoxmutt $hoxstep $hoxhnmut
       ./run.e $name $ini $mut $gen $i $mutstep $pleio $hox $hoxmutt $hoxstep $hoxhnmut #& no "&" ... so this is NOT happening in paralel ... #making the mutant files / no psize needed here, size of file is the psize
   done
  echo "done mutating" 
  #Safe genealogy
  cat *"genealoo"* >> genealogy.dat #copy all teeth to genealogy, even if no mutation
  cat *"TORUN1" > TORUN.txt      #all go TORUN, then check if mutated (A = it has mutated, B= it has not mutated)
  cat *"genealoo.t00t" > nonmuta.tmp #no mutation in either tooth, writen in evol.f90
  #put in fitness_now.dat the fitness of the non mutated ones, fitness is the same than their parent
  #the file nonmuta.tmp contains the non-mutated files that are included in current generation      
  x=$(wc -l nonmuta.tmp | awk '{ print $1 }')
  echo $x
  for i in $(seq 1 $x); do
   #   echo $gen
        name1=$(awk -F" " 'FNR == '$i' {print $1}' nonmuta.tmp) #name of parental 1
        name2=$(awk -F" " 'FNR == '$i' {print $2}' nonmuta.tmp) #name of parental 2
        name3=$(awk -F" " 'FNR == '$i' {print $3}' nonmuta.tmp) #name of offspring 1
        name4=$(awk -F" " 'FNR == '$i' {print $4}' nonmuta.tmp) #name of offspring 2
        fit=$(awk /$name1/'{ print $5; exit}' fitness_all.dat)  #fitness of parental, always in pairs, finding "half" a parental (one teeth) is enough to find fitness
        echo ${name1}profile__up.txt $name3 ${name2}profile_down.txt $name4  $fit $gen >> fitness_now.dat 
  done     
  #**********************************************************
  
  rm  *"genealoo"* nonmuta.tmp
  rm *"TORUN1"
  ini=1 #now changing the files with mutations, not just introducing some initial noise in parameters
  
  psize2=$(wc -l TORUN.txt | awk '{ print $1 }') #number of teeth to run (that have mutated)  !teeth in same mouth do not run in paralel right now TODO?
  i=0

  #Develop the tooth
  while [ $i -lt $psize2 ]
   do
     cpuRUNNING=$(ps -aef | grep runt$1.e | wc -l)
     if [ "$cpuRUNNING" -le "$maxCPU2" ] ; then
       i=`expr $i + $one`
       name1=$(awk -F" " 'FNR == '$i' {print $1}' TORUN.txt) #name of parental 1
       name2=$(awk -F" " 'FNR == '$i' {print $2}' TORUN.txt) #name of parental 2
       mut1a=$(awk -F" " 'FNR == '$i' {print $3}' TORUN.txt) #name of offspring 1
       mut1b=$(awk -F" " 'FNR == '$i' {print $4}' TORUN.txt) #has offspring 1 mutated? A yes, B no
       mut2a=$(awk -F" " 'FNR == '$i' {print $5}' TORUN.txt) #name of offspring 2
       mut2b=$(awk -F" " 'FNR == '$i' {print $6}' TORUN.txt) #has offspring 2 mutated? A yes, B no
       echo $i "running" $name1 $name2 $mut1a $mut1b $mut2a $mut2b
       #timeout -k 1 5m ./runt$1.e $name1 $name2 $mut1a $mut1b $mut2a $mut2b >out$i 2> out2_$i & #running development      fi
       timeout -k 1 3m ./runt$1.e $name1 $name2 $mut1a $mut1b $mut2a $mut2b >out$i 2> out2_$i || echo "Development_too_long" > $name1 $name2 > $name1$name2".toolong" & #running development 
     fi
   done
   
    check=0  #keep checking until all done
    while [ $check -lt $one ]
     do     
       cpuRUNNING=$(ps -aef | grep runt$1.e | wc -l)
       if [ "$cpuRUNNING" -le "$one" ] ; then # development of all teeth has concluded      
         check=1  #allows to leave loop
         i=0    
         # make list of teeth to meassure fitness using files writen at the end of teeth devo (in source_original.f90)
          cat *_info.txt > TOFIT.txt             
         psize2=$(wc -l TOFIT.txt | awk '{ print $1 }') #number of teeth to calculate fitness
         #calculate fitness      # TODO modificar occlusion para q lea los inputs  / hacer lista para la siguiente gen, esta lista debe contener n=psize segun fitness
         while [ $i -lt $psize2 ]
         do
          cpuRUNNING=$(ps -aef | grep fitness$1.e | wc -l)
    	   if [ "$cpuRUNNING" -le "$maxCPU" ] ; then
    	     i=`expr $i + $one`
             name1=$(awk -F" " 'FNR == '$i' {print $1}' TOFIT.txt)
             name2=$(awk -F" " 'FNR == '$i' {print $2}' TOFIT.txt)
             name3=$(awk -F" " 'FNR == '$i' {print $3}' TOFIT.txt)
             name4=$(awk -F" " 'FNR == '$i' {print $4}' TOFIT.txt)
	     echo $i "fit"
             ./fitness$1.e $name1 $name2 $name3 $name4 $gen $fittype2 >outfit$i 2> out2fit_$i & #running occlusion 
           fi
         done                
       fi
     done  
     #wait for fitness calculation to be done
     check=0 #keep checking until all done
     while [ $check -lt $one ]
     do     
       cpuRUNNING=$(ps -aef | grep fitness$1.e | wc -l)
       #echo "hola!!!!!!!!" "$cpuRUNNING"
       if [ "$cpuRUNNING" -le "$one" ] ; then           
           cat *FITNESS* >> fitness_now.dat  #only for current generation        
           cat fitness_now.dat >> fitness_all.dat #safe data all generations                            
           rm *FITNESS*
           rm *_info.txt
           check=1     
       fi
     done 
   
   #choose the one that will contribute to the next generation
   #50%
    
    sort -gk 5,5 fitness_now.dat > fitness_now.tmp #order file using the 5th column, which contains the distance values, lower distance is higher fitness 
    head -$half fitness_now.tmp > top50p.dat
#exit
    #making sure at that psize stays the same. Due to many abortions it can happen that not enough individuals. In that case we take more starting with the ones with better fitness
 
   n=$(wc -l top50p.dat | awk '{ print $1 }')
   while [ $n -lt $half ]
   do      
     n=$(wc -l top50p.dat | awk '{ print $1 }')
     mis=`expr $half - $n` 
     head -$mis fitness_now.tmp >> top50p.dat		      
   done
   n=$(wc -l top50p.dat | awk '{ print $1 }')
   if [ $n -gt $half ] ; then
         head -$half top50p.dat > top50p.dat
   fi

   echo | awk -F" " '{print $2, $4}' top50p.dat > TOMUT.txt
   echo | awk -F" " '{print $2, $4}' top50p.dat >> TOMUT.txt #duplicate, since we chose the top 50%, every individual will have two offsprings
   
   #safe previous generation individuals to calculate fitness for when change in selection happens
   cp fitness_now.dat fitness_1gen_before.dat 
   
   
   rm fitness_now.tmp fitness_now.dat   
   
   
   mkdir $gen
   mv *.off $gen
   mv *.xyz $gen
   tar -czf $gen.tar.gz $gen
   rm -r $gen
   
done 

echo "Done!" 

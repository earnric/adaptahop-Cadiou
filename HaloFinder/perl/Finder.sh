for i in {1..51}
do
     sed -i s/HaloMaker/HaloMaker_$i/g Halo.job
     sed -i s/j=PBS_ARRAYID/j=$i/g Halo.job
     qsub Halo.job
     sed -i s/HaloMaker_$i/HaloMaker/g Halo.job
     sed -i s/j=$i/j=PBS_ARRAYID/g Halo.job
done
 

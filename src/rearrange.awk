#! /bin/awk -f

#expects step, ncontrols, ncohort

BEGIN{
    chr=1;
    pos=2;
    id=3
    ref=4;
    alt=5;
    svlen=6;
    svtype=7;
    collapseid=8;
    svgt=9;
    svvarreads=10;
    svrefreads=11
    totalreads=12 #at end
    qend=11
    if(ncohort>1){
        qstart=9
        qend=qstart+(ncohort*3)-1
    }
    controlstart=qend+1
    controlend=controlstart+(ncontrols*3)-1
    genotypes=13
    OFS="\t"
}
{
    if(ncohort==1){
        tr=$svvarreads+$svrefreads
        $svrefreads=$svrefreads"\t"tr
    }else{
        for (i=qstart;i<=(qend-2);i+=step){
            tr=$(i+1)+$(i+2)
            $(i+2)=$(i+2)"\t"tr
        }
    }
    for(i=controlstart;i<=(controlend-2);i+=step){
        $(i+1)=""
        $(i+2)=""
    }
    $2=$2"\t."
    print $0
    #printf("%s\t%s\t%s\t%s\t%s\t%s",
    #     $chr, $pos, $ref, $alt, $svlen, $svtype, $svgt)
    # Now print the variable columns in a loop
    #for (i = qstart; i <= qend; i += step) {
    #    totalreads=($i+1)+($i+2)
    #    printf("\t%s\t%s\t%s\t%s", $i, $i+1, $i+2, totalreads)
    #}
    # And print the control variable columns in a loop (gt only)
    #for (i = start; i <= end; i += step) {
    #    printf("\t%s", $i)
    #}
    #printf("\n")
}

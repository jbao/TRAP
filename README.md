TRAP
====

Modified code from the TF affinity prediction

See also http://www.molgen.mpg.de/~manke/papers/TFaffinities/

Batch processing
----------------

In order to analyze multiple sequences, one has to first extract the promoter 
sequences, eg.
    
    human.promoter <- read.DNAStringSet('~/data/fimo/human_upstream1000.fa')

where `human_upstream1000.fa` is the whole genome promoter sequence file, one 
can then pick her genes of interest (eg. all upregulated genes) and find out
the ENTREZ IDs
    
    tf.entrez <- unlist(mget(as.character(df.uptime$gene), org.Hs.egSYMBOL2EG, 
        ifnotfound=NA))
    tf.entrez <- tf.entrez[!is.na(tf.entrez)]

which have to be converted to REFSEQ IDs to match the promoter file

    tf.refseq <- unlist(mget(tf.entrez, org.Hs.egREFSEQ))
    tf.refseq <- tf.refseq[grep('NM',tf.refseq)]
    idx <- grep(paste(tf.refseq,collapse='_|'), names(human.promoter))

Here we define a function to convert REFSEQ IDs back to gene symbols, such 
that our sequence files will have the consistent gene names

    refseq2symbol <- function(refseq.str) {
        all.str <- unlist(strsplit(refseq.str, '_'))
        refseq <- paste(all.str[1:2], collapse='_')
        entrez <- unlist(mget(refseq,org.Hs.egREFSEQ2EG,ifnotfound=NA))
        if (is.na(entrez))
            return('')
        else
            return(get(entrez, org.Hs.egSYMBOL))
    }

Now we can define our output directory and dump the sequence files into it

    wd <- '~/data/hdf/qpcr_seq/'
    for (i in 1:length(idx)) {
        symbol <- refseq2symbol(names(human.promoter)[idx[i]])
        if (symbol != '') {
            seq.file <- paste(wd,symbol,'.fa',sep='')
            write.XStringSet(human.promoter[idx[i]], seq.file)
        }
        symbol
    }

With all that in place, we can then execute the TRAP program on all sequence 
files (eg. using the parallel functionality in `foreach`)

    all.file <- dir(wd, 'fa')
    foreach(i=1:length(all.file)) %dopar% {
        seq.file <- paste(wd,all.file[i],sep='')
        res.file <- paste(wd,sub('fa','trap',all.file[i]),sep='')
        system(paste('./TRAP ~/data/TRANSFAC/TFP_2012.3/dat/matrix.dat',seq.file,
            '>',res.file))
    }

where the `matrix.dat` is the position score matrix provided by TRANSFAC, and 
the result files will have the extension of `.trap`, while locating in the same
directory. The result files can then be individually processed to calculate the
normalized affinity score.

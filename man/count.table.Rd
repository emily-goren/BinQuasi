% File src/library/QuasiSeq/man/count.table.Rd
\name{count.table}
\alias{count.table}
\title{Create a Matrix of ChIP-seq Count Data}
\description{Create a matrix of ChIP-seq count data from sorted bam
  files using a non-overlapping genomic partition.}
\usage{
count.table(dir, ChIP.files, control.files, bin.size = NULL,
           frag.length = NULL, minimum.count = 20)
}
\arguments{
  \item{dir}{Directory where the sorted bam files are saved.}
  \item{ChIP.files}{File names (with file extensions) of the ChIP sample
  files in sorted bam format.}
  \item{control.files}{File names (with file extensions) of the
  input/control sample files in sorted bam format.}
  \item{bin.size}{Window size, constant across all samples, used to
    generate a non-overlapping partition for
    counts. If \code{NULL}, an estimate will be used (see details).}
  \item{frag.length}{Average length of the ChIP fragments in each sample
    provided. Reads are extended to this length from their 3' ends. If \code{NULL}, cross correlation will be used to
    estimate the fragment length of each sample (see details).}
  \item{minimum.count}{The count threshold used for filtering out
    windows with sparse counts. Any genomic window with counts less than
  this value accross all samples will be removed.}
}
\value{A list containing:
  \item{counts}{data.frame with rows corresponding to genomic windows and
       columns for the chromosomes, start and end locations, as well as
       a column for the counts of each sample.}
     \item{bin.size}{the bin.size used to create the counts dataframe.}
     \item{fragment.length}{vector of the fragment lengths used to
       extend the reads in each sample.}
     \item{filter}{the count threshold used to create the counts
       dataframe. Windows with counts summed accross all samples that fall below this value were removed.}
     }
\details{
This function creates a count table of ChIP sequencing data
(supplied as sorted bam files) using a non-overlapping partition accross
the genome.

The fragment length (if not provided) is estimated using the
cross-correlation method of Ramachandran et al. (2013). A fregment
length is estimated for each sample, after removing dupliate reads, by taking the average over
all chromosomes in the sample. Estimation is performed at 5 bp
resolution and restricted to a minimum fragment length of 50 bp and maximum of 600 bp.

The bin size (if not provided) is selected using a procedure by
Shimazaki and Shinomoto (2007) based on minimzing the Mean-Integrated
Squared Error for a time-dependent Poisson point process. This procedure
is applied to each ChIP sample (at 5 bp resolution, restricted to a
minimum of 50 bp and maximum of 1000 bp), and the minimum across all ChIP samples
is returned as the bin size.

For a given sample and window, the count is determined as the number of
fragments overlapping the window.
}

\references{
Shimazaki and Shinomoto (2007)  "A method for selecting the bin size of a
time histogram" \emph{Neural computation}, \bold{19}(6), 1503-27.

Ramachandran, Palidwor, Porter,  and Perkins (2013)
 "MaSC: mappability-sensitive cross-correlation for estimating mean fragment length of single-end short-read sequencing data" \emph{Bioinformatics} \bold{29}(4), 444-50.
}

\author{Emily Goren \email{emily.goren@gmail.com}}
\examples{
\dontrun{

    d <- count.table(dir = '/Users/username/mybamfolder/',
                     ChIP.files = c('chip1.bam', 'chip2.bam'),
                     control.files = c('inp1.bam', 'inp2.bam'),
                     bin.size = 60, frag.length = c(101, 300, 150, 10),
                     minimum.count = 20)
    head(d$counts)
}
}
\keyword{ChIP-seq, quasi-likelihood, next-generation sequencing}


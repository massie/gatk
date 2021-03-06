package org.broadinstitute.hellbender.engine.spark.datasources;

import com.google.api.services.storage.Storage;
import com.google.cloud.genomics.dataflow.readers.bam.BAMIO;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.mapreduce.Job;
import org.apache.parquet.avro.AvroParquetInputFormat;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.bdgenomics.formats.avro.AlignmentRecord;
import org.broadinstitute.hellbender.engine.AuthHolder;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.BDGAlignmentRecordToGATKReadAdapter;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.SAMRecordToGATKReadAdapter;
import org.seqdoop.hadoop_bam.AnySAMInputFormat;
import org.seqdoop.hadoop_bam.SAMRecordWritable;
import org.seqdoop.hadoop_bam.util.SAMHeaderReader;

import java.io.IOException;
import java.io.Serializable;
import java.util.List;

/** Loads the reads from disk either serially (using samReaderFactory) or in parallel using Hadoop-BAM.
 * The parallel code is a modified version of the example writing code from Hadoop-BAM.
 */
public class ReadsSparkSource implements Serializable {
    private static final long serialVersionUID = 1L;
    public static final int DEFAULT_SPLIT_SIZE = 10485760;

    private transient final JavaSparkContext ctx;
    public ReadsSparkSource(JavaSparkContext ctx) {
        this.ctx = ctx;
    }


    /**
     * Loads Reads using Hadoop-BAM. For local files, bam must have the fully-qualified path,
     * i.e., file:///path/to/bam.bam.
     * @param bam file to load
     * @param intervals intervals of reads to include.
     * @return RDD of (SAMRecord-backed) GATKReads from the file.
     */
    public JavaRDD<GATKRead> getParallelReads(final String bam, final List<SimpleInterval> intervals) {
        return getParallelReads(bam, intervals, DEFAULT_SPLIT_SIZE);
    }

    /**
     * Loads Reads using Hadoop-BAM. For local files, bam must have the fully-qualified path,
     * i.e., file:///path/to/bam.bam.
     * @param bam file to load
     * @param intervals intervals of reads to include.
     * @param splitSize maximum bytes of bam file to read into a single partition, increasing this will result in fewer partitions
     * @return RDD of (SAMRecord-backed) GATKReads from the file.
     */
    public JavaRDD<GATKRead> getParallelReads(final String bam, final List<SimpleInterval> intervals, long splitSize) {
        Configuration conf = new Configuration();
        // reads take more space in memory than on disk so we need to limit the split size
        conf.set("mapreduce.input.fileinputformat.split.maxsize", Long.toString(splitSize));

        JavaPairRDD<LongWritable, SAMRecordWritable> rdd2 = ctx.newAPIHadoopFile(
                bam, AnySAMInputFormat.class, LongWritable.class, SAMRecordWritable.class,
                conf);

        return rdd2.map(v1 -> {
            SAMRecord sam = v1._2().get();
            if (samRecordOverlaps(sam, intervals)) {
                try {
                    return (GATKRead) new SAMRecordToGATKReadAdapter(sam);
                } catch (SAMException e) {
                    // TODO: add stringency
                }            
            }
            return null;

        }).filter(v1 -> v1 != null);
    }

    /**
     * Loads Reads using Hadoop-BAM. For local files, bam must have the fully-qualified path,
     * i.e., file:///path/to/bam.bam. This excludes unmapped reads.
     * @param bam file to load
     * @return RDD of (SAMRecord-backed) GATKReads from the file.
     */
    public JavaRDD<GATKRead> getParallelReads(final String bam) {
        return getParallelReads(bam, DEFAULT_SPLIT_SIZE);
    }

    /**
     * Loads Reads using Hadoop-BAM. For local files, bam must have the fully-qualified path,
     * i.e., file:///path/to/bam.bam. This excludes unmapped reads.
     * @param bam file to load
     * @param splitSize maximum bytes of bam file to read into a single partition, increasing this will result in fewer partitions
     * @return RDD of (SAMRecord-backed) GATKReads from the file.
     */
    public JavaRDD<GATKRead> getParallelReads(final String bam, int splitSize) {
        final SAMFileHeader readsHeader = getHeader(ctx, bam, null);
        List<SimpleInterval> intervals = IntervalUtils.getAllIntervalsForReference(readsHeader.getSequenceDictionary());
        return getParallelReads(bam, intervals, splitSize);
    }

    /**
     * Loads ADAM reads stored as Parquet.
     * @param inputPath path to the Parquet data
     * @return RDD of (ADAM-backed) GATKReads from the file.
     */
    public JavaRDD<GATKRead> getADAMReads(final String inputPath, final List<SimpleInterval> intervals, final SAMFileHeader header) throws IOException {
        Job job = Job.getInstance(ctx.hadoopConfiguration());
        AvroParquetInputFormat.setAvroReadSchema(job, AlignmentRecord.getClassSchema());
        Broadcast<SAMFileHeader> bHeader;
        if (header == null) {
            bHeader= ctx.broadcast(null);
        } else {
            bHeader = ctx.broadcast(header);
        }
        @SuppressWarnings("unchecked")
        JavaRDD<AlignmentRecord> recordsRdd = ctx.newAPIHadoopFile(
                inputPath, AvroParquetInputFormat.class, Void.class, AlignmentRecord.class, job.getConfiguration())
                .values();
        JavaRDD<GATKRead> readsRdd = recordsRdd.map(record -> new BDGAlignmentRecordToGATKReadAdapter(record, bHeader.getValue()));
        JavaRDD<GATKRead> filteredRdd = readsRdd.filter(record -> samRecordOverlaps(record.convertToSAMRecord(header), intervals));
        return filteredRdd;
    }

    /**
     * Loads the header using Hadoop-BAM.
     * @param filePath path to the bam.
     * @param auth authentication information if using GCS.
     * @return the header for the bam.
     */
    public static SAMFileHeader getHeader(final JavaSparkContext ctx, final String filePath, final AuthHolder auth) {
        final SAMFileHeader samFileHeader;

        // GCS case
        if (BucketUtils.isCloudStorageUrl(filePath)) {
            try {
                Storage.Objects storageClient = auth.makeStorageClient();
                try (final SamReader reader = BAMIO.openBAM(storageClient, filePath, ValidationStringency.DEFAULT_STRINGENCY)) {
                    return reader.getFileHeader();
                }
            } catch (Exception e) {
                throw new GATKException("Failed to read bams header from " + filePath + ".", e);
            }
        }

        try {
             samFileHeader = SAMHeaderReader.readSAMHeaderFrom(new Path(filePath), ctx.hadoopConfiguration());
        } catch (IOException e) {
            throw new GATKException("unable to loader header: " + e);
        }
        return samFileHeader;
    }

    /**
     * Tests if a given SAMRecord overlaps any interval in a collection.
     */
    //TODO: remove this method when https://github.com/broadinstitute/hellbender/issues/559 is fixed
    private static boolean samRecordOverlaps(final SAMRecord record, final List<SimpleInterval> intervals ) {
        if (intervals == null || intervals.isEmpty()) {
            return true;
        }
        for (SimpleInterval interval : intervals) {
            if (record.getReadUnmappedFlag() && record.getAlignmentStart() != SAMRecord.NO_ALIGNMENT_START) {
                // This follows the behavior of htsjdk's SamReader which states that "an unmapped read will be returned
                // by this call if it has a coordinate for the purpose of sorting that is in the query region".
                int start = record.getAlignmentStart();
                return interval.getStart() <= start && interval.getEnd() >= start;
            } else  if (interval.overlaps(record)) {
                return true;
            }
        }
        return false;
    }
}

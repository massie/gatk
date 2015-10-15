package org.broadinstitute.hellbender.tools.spark.pipelines;

import com.google.cloud.dataflow.sdk.transforms.SerializableFunction;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.ReadContextData;
import org.broadinstitute.hellbender.engine.spark.AddContextDataToReadSpark;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.JoinStrategy;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.engine.spark.datasources.VariantsSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.ApplyBQSRArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.transforms.ApplyBQSRSparkFn;
import org.broadinstitute.hellbender.tools.spark.transforms.BaseRecalibratorSparkFn;
import org.broadinstitute.hellbender.tools.walkers.bqsr.BaseRecalibrator;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;
import org.broadinstitute.hellbender.utils.recalibration.BaseRecalibrationEngine;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationArgumentCollection;
import org.broadinstitute.hellbender.utils.recalibration.RecalibrationReport;
import org.broadinstitute.hellbender.utils.variant.Variant;

import java.io.IOException;
import java.util.List;
@CommandLineProgramProperties(
        summary = "This tools performs 2 steps of BQSR - crection of recalibration tables and rewriting of the bam, without writing the tables to disk. ",
        oneLineSummary = "Both steps of BQSR",
        usageExample = "BQSRSpark -I in.bam --BQSRKnownVariants in.vcf -O out.bam",
        programGroup = SparkProgramGroup.class
)
/**
 * BQSR. The final result is analysis-ready reads.
 */
public final class BQSRSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Override
    public boolean requiresReads() { return true; }

    @Override
    public boolean requiresReference() { return true; }

    @Argument(doc = "the known variants", shortName = "BQSRKnownVariants", fullName = "baseRecalibrationKnownVariants", optional = true)
    protected List<String> baseRecalibrationKnownVariants;

    @Argument(doc = "the output bam", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    protected String output;

    @Argument(doc = "If specified, shard the output bam", shortName = "shardedOutput", fullName = "shardedOutput", optional = true)
    private boolean shardedOutput = false;

    @Argument(doc = "the join strategy for reference bases", shortName = "joinStrategy", fullName = "joinStrategy", optional = true)
    private JoinStrategy joinStrategy = JoinStrategy.SHUFFLE;
    
    @Override
    public SerializableFunction<GATKRead, SimpleInterval> getReferenceWindowFunction() {
        return BaseRecalibrationEngine.BQSR_REFERENCE_WINDOW_FUNCTION;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        JavaRDD<GATKRead> markedReads = getReads();
        VariantsSparkSource variantsSparkSource = new VariantsSparkSource(ctx);

        if ( baseRecalibrationKnownVariants.isEmpty() ) { // Warn the user if no dbSNP file or other variant mask was specified
            throw new UserException.CommandLineException(BaseRecalibrator.NO_DBSNP_EXCEPTION);
        }

        // TODO: workaround for known bug in List version of getParallelVariants
        if ( baseRecalibrationKnownVariants.size() > 1 ) {
            throw new GATKException("Cannot currently handle more than one known sites file, " +
                                    "as getParallelVariants(List) is broken");
        }
        JavaRDD<Variant> bqsrKnownVariants = variantsSparkSource.getParallelVariants(baseRecalibrationKnownVariants.get(0));

        // TODO: Look into broadcasting the reference to all of the workers. This would make AddContextDataToReadSpark
        // TODO: and ApplyBQSRStub simpler (#855).
        JavaPairRDD<GATKRead, ReadContextData> rddReadContext = AddContextDataToReadSpark.add(markedReads, getReference(), bqsrKnownVariants, joinStrategy);
        // TODO: broadcast the reads header?

        //note: we use the reference dictionary from the reads themselves.
        final RecalibrationReport bqsrReport = BaseRecalibratorSparkFn.apply(rddReadContext, getHeaderForReads(), getHeaderForReads().getSequenceDictionary(), new RecalibrationArgumentCollection());
        final Broadcast<RecalibrationReport> reportBroadcast = ctx.broadcast(bqsrReport);
        final JavaRDD<GATKRead> finalReads = ApplyBQSRSparkFn.apply(markedReads, reportBroadcast, getHeaderForReads(), new ApplyBQSRArgumentCollection());

        try {
            ReadsSparkSink.writeReads(ctx, output, finalReads, getHeaderForReads(), shardedOutput ? ReadsWriteFormat.SHARDED : ReadsWriteFormat.SINGLE);
        } catch (IOException e) {
            throw new GATKException("unable to write bam: " + e);
        }
    }
}

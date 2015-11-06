package org.broadinstitute.hellbender.tools.spark.pipelines;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.Arrays;

public class BQSRPipelineSparkIntegrationTest extends CommandLineProgramTest {
    private static final class BQSRTest {
        final String referenceURL;
        final String bam;
        final String knownSites;
        final String args;
        final String expectedFileName;

        private BQSRTest(String referenceURL, String bam, String knownSites, String args, String expectedFileName) {
            this.referenceURL = referenceURL;
            this.bam = bam;
            this.knownSites = knownSites;
            this.args = args;
            this.expectedFileName = expectedFileName;
        }

        public String getCommandLineNoApiKey() {
            return  " -R " + referenceURL +
                    " -I " + bam +
                    " " + args +
                    (knownSites.isEmpty() ? "": " --baseRecalibrationKnownVariants " + knownSites) +
                    " -O %s";
        }

        public String getCommandLine() {
            return  getCommandLineNoApiKey() +
                    " --apiKey " + getGCPTestApiKey();
        }

        @Override
        public String toString() {
            return String.format("BQSR(bam='%s', args='%s')", bam, args);
        }
    }

    private String getResourceDir(){
        return getTestDataDir() + "/" + "BQSR" + "/";
    }

    @DataProvider(name = "BQSRLocalRefTest")
    public Object[][] createBQSRLocalRefTestData() {
        final String GRCh37Ref_2021 = b37_reference_20_21;
        final String GRCh37Ref2bit_chr2021 = b37_2bit_reference_20_21;
        final String hiSeqBam_chr20 = getResourceDir() + "CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.bam";
        final String dbSNPb37_2021 = getResourceDir() + "dbsnp_138.b37.excluding_sites_after_129.ch20.1m-1m1k.vcf";

        return new Object[][]{
                // input local, computation local.
                {new BQSRTest(GRCh37Ref2bit_chr2021, hiSeqBam_chr20, dbSNPb37_2021, "--joinStrategy BROADCAST", getResourceDir() + "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.recalibrated.bam")},
                {new BQSRTest(GRCh37Ref_2021, hiSeqBam_chr20, dbSNPb37_2021, "--joinStrategy SHUFFLE", getResourceDir() + "expected.CEUTrio.HiSeq.WGS.b37.ch20.1m-1m1k.NA12878.recalibrated.bam")},
        };
    }

    @Test(dataProvider = "BQSRLocalRefTest")
    public void testBQSRLocalRef(BQSRTest params) throws IOException {
        ArgumentsBuilder ab = new ArgumentsBuilder().add(params.getCommandLineNoApiKey());
        IntegrationTestSpec spec = new IntegrationTestSpec(
                ab.getString(),
                Arrays.asList(params.expectedFileName));
        spec.executeTest("testBQSR-" + params.args, this);
    }
}

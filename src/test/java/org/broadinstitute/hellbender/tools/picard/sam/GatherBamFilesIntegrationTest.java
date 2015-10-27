package org.broadinstitute.hellbender.tools.picard.sam;

import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.cram.build.CramIO;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public final class GatherBamFilesIntegrationTest extends CommandLineProgramTest {
    private static final File TEST_DATA_DIR = new File(getTestDataDir(), "picard/sam/GatherBamFiles");
    private static final File ORIG_BAM = new File(TEST_DATA_DIR, "orig.bam");
    private static final List<File> SPLIT_BAMS = Arrays.asList(
            new File(TEST_DATA_DIR, "indUnknownChrom.bam"),
            new File(TEST_DATA_DIR, "indchr1.bam"),
            new File(TEST_DATA_DIR, "indchr2.bam"),
            new File(TEST_DATA_DIR, "indchr3.bam"),
            new File(TEST_DATA_DIR, "indchr4.bam"),
            new File(TEST_DATA_DIR, "indchr5.bam"),
            new File(TEST_DATA_DIR, "indchr6.bam"),
            new File(TEST_DATA_DIR, "indchr7.bam"),
            new File(TEST_DATA_DIR, "indchr8.bam")
    );

    public String getTestedClassName() {
        return GatherBamFiles.class.getSimpleName();
    }

    private void testTheGathering(final File aggregateFile, final List<File> splitFiles, final File referenceFile, final String outputExtension) throws Exception {
        final File outputFile = BaseTest.createTempFile("gatherBamFilesTest.samFile.", outputExtension);
        final List<String> args = new ArrayList<>();
        for (final File splitBam : splitFiles) {
            args.add("--INPUT");
            args.add(splitBam.getAbsolutePath());
        }
        if (null != referenceFile) {
            args.add("--R");
            args.add(referenceFile.getAbsolutePath());
        }
        args.add("--OUTPUT");
        args.add(outputFile.getAbsolutePath());
        runCommandLine(args);
        SamAssertionUtils.assertSamsEqual(aggregateFile, outputFile, referenceFile);
        SamAssertionUtils.assertSamsNonEqual(aggregateFile, splitFiles.get(0), referenceFile); // sanity check
    }

    @Test
    public void testTheBAMGathering() throws Exception {
        testTheGathering(ORIG_BAM, SPLIT_BAMS, null, BamFileIoUtils.BAM_FILE_EXTENSION);
    }

    private static final File ORIG_CRAM = new File(TEST_DATA_DIR, "orig.cram");
    private static final List<File> SPLIT_CRAMS = Arrays.asList(
            // these need to be presented in coordinate order with unknown last
            new File(TEST_DATA_DIR, "indchr1.cram"),
            new File(TEST_DATA_DIR, "indchr2.cram"),
            new File(TEST_DATA_DIR, "indchr3.cram"),
            new File(TEST_DATA_DIR, "indchr4.cram"),
            new File(TEST_DATA_DIR, "indchr5.cram"),
            new File(TEST_DATA_DIR, "indchr6.cram"),
            new File(TEST_DATA_DIR, "indchr7.cram"),
            new File(TEST_DATA_DIR, "indchr8.cram"),
            new File(TEST_DATA_DIR, "indUnknownChrom.cram")
    );

    @Test
    public void testTheCRAMGathering() throws Exception {
        testTheGathering(ORIG_CRAM, SPLIT_CRAMS, new File(TEST_DATA_DIR, "basic.fasta"), CramIO.CRAM_FILE_EXTENSION);
    }

    private static final List<File> INVALID_SPLIT_CRAMS = Arrays.asList(
            // present the files out of order
            new File(TEST_DATA_DIR, "indUnknownChrom.cram"),
            new File(TEST_DATA_DIR, "indchr1.cram")
    );

    @Test(expectedExceptions=IllegalArgumentException.class)
    public void testTheInvalidSortCRAMGathering() throws Exception {
        testTheGathering(ORIG_CRAM, INVALID_SPLIT_CRAMS, new File(TEST_DATA_DIR, "basic.fasta"), ".cram");
    }
}

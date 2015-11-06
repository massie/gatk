package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;

import java.io.IOException;
import java.util.*;

@CommandLineProgramProperties(summary="Gather clustered split reads using spark",
        oneLineSummary="Gather clustered split reads using spark",
        programGroup = SparkProgramGroup.class)
public class GatherSplitReadsSpark extends GATKSparkTool
{
    private static final long serialVersionUID = 1l;

    private static final int MIN_SOFT_CLIP_LEN = 30; // minimum length of an interesting soft clip
    private static final byte MIN_QUALITY = 15; // minimum acceptable quality in a soft-clip window
    private static final int MAX_LOCUS_DIST = 500; // stale evidence distance, ought to be somewhat longer than a read
    private static final int CLUSTER_WINDOW = 2; // size of locus window in which we cluster event evidence
    private static final int MIN_EVIDENCE = 3; // minimum evidence count in a cluster

    @Argument(doc = "the output bam", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String output;

    @Override
    public boolean requiresReads()
    {
        return true;
    }

    private static final class EventEvidence implements Comparable<EventEvidence>
    {
        public EventEvidence( final int someContigIdx, final int someLocus )
        {
            contigIdx = someContigIdx;
            locus = someLocus;
            name = "";
        }

        public EventEvidence( final int someContigIdx, final int someLocus, final String someName )
        {
            contigIdx = someContigIdx;
            locus = someLocus;
            name = someName;
        }

        public int getContigIdx() { return contigIdx; }
        public int getLocus() { return locus; }

        public int compareTo( final EventEvidence ee )
        {
            int result = Integer.compare(contigIdx,ee.contigIdx);
            if ( result == 0 ) result = Integer.compare(locus,ee.locus);
            if ( result == 0 ) result = name.compareTo(ee.name);
            return result;
        }

        private final int contigIdx;
        private final int locus;
        private final String name;
    }

    private static final class SplitReadClusterer implements Iterator<GATKRead>, Iterable<GATKRead>
    {
        public SplitReadClusterer( Iterator<GATKRead> someIterator, SAMSequenceDictionary dictionary )
        {
            this.inputIterator = someIterator;
            this.dictionary = dictionary;
        }

        @Override
        public boolean hasNext()
        {
            while ( !outputIterator.hasNext() )
            {
                if ( !inputIterator.hasNext() )
                    return false;
                outputIterator = processRead(inputIterator.next());
            }
            return true;
        }

        @Override
        public GATKRead next()
        {
            if ( !hasNext() )
                throw new NoSuchElementException("Iterator<GATKRead> is exhausted.");
            return outputIterator.next();
        }

        @Override
        public Iterator<GATKRead> iterator() { return this; }

        private Iterator<GATKRead> processRead( GATKRead read )
        {
            if (!read.failsVendorQualityCheck() && !read.isUnmapped() && read.getMappingQuality() > 0) {
                final List<CigarElement> cigarElements = read.getCigar().getCigarElements();
                final ListIterator<CigarElement> itr0 = cigarElements.listIterator();
                if (itr0.hasNext()) {
                    CigarElement firstEle = itr0.next();
                    if (firstEle.getOperator() == CigarOperator.HARD_CLIP && itr0.hasNext())
                        firstEle = itr0.next();
                    if (firstEle.getOperator() == CigarOperator.SOFT_CLIP &&
                            firstEle.getLength() >= MIN_SOFT_CLIP_LEN &&
                            highQualityRegion(read.getBaseQualities(), 0))
                        return cluster(read, read.getStart()); // NON-STRUCTURED return
                }
                final ListIterator<CigarElement> itrN = cigarElements.listIterator(cigarElements.size());
                if (itrN.hasPrevious()) {
                    CigarElement lastEle = itrN.previous();
                    if (lastEle.getOperator() == CigarOperator.HARD_CLIP && itrN.hasPrevious())
                        lastEle = itrN.previous();
                    if (lastEle.getOperator() == CigarOperator.SOFT_CLIP &&
                            lastEle.getLength() >= MIN_SOFT_CLIP_LEN &&
                            highQualityRegion(read.getBaseQualities(), read.getLength() - lastEle.getLength()))
                        return cluster(read, read.getEnd()); // NON-STRUCTURED return
                }
            }
            return Collections.emptyIterator();
        }

        private Iterator<GATKRead> cluster( final GATKRead read, final int locus )
        {
            final int contigIdx = dictionary.getSequenceIndex(read.getContig());
            final EventEvidence newEvidence = new EventEvidence(contigIdx,locus,read.getName());
            locMap.put(newEvidence,read);

            // clean out old stuff that can't possibly be interesting anymore
            final Iterator<Map.Entry<EventEvidence,GATKRead>> itr = locMap.entrySet().iterator();
            while ( itr.hasNext() )
            {
                final EventEvidence oldEvidence = itr.next().getKey();
                if ( oldEvidence.getContigIdx() == newEvidence.getContigIdx() &&
                        oldEvidence.getLocus() + MAX_LOCUS_DIST > newEvidence.getLocus() )
                    break;
                itr.remove();
            }

            // find all the evidence in a window surrounding the locus of the new evidence
            final SortedMap<EventEvidence,GATKRead> windowMap = locMap
                    .tailMap(new EventEvidence(newEvidence.getContigIdx(), newEvidence.getLocus() - CLUSTER_WINDOW))
                    .headMap(new EventEvidence(newEvidence.getContigIdx(), newEvidence.getLocus() + CLUSTER_WINDOW + 1));
            if ( windowMap.size() >= MIN_EVIDENCE )
            {
                final List<GATKRead> list = new LinkedList<>();
                for ( Map.Entry<EventEvidence, GATKRead> entry : windowMap.entrySet() )
                {
                    if (entry.getValue() != null)
                    {
                        list.add(entry.getValue());
                        entry.setValue(null);
                    }
                }
                return list.iterator(); // NON-STRUCTURED return
            }
            return Collections.emptyIterator();
        }

        private static boolean highQualityRegion( final byte[] quals, int idx )
        {
            for ( final int end = idx+MIN_SOFT_CLIP_LEN; idx != end; ++idx )
                if ( quals[idx] < MIN_QUALITY )
                    return false; // NON-STRUCTURED return
            return true;
        }

        private final Iterator<GATKRead> inputIterator;
        private final SAMSequenceDictionary dictionary;
        private final SortedMap<EventEvidence,GATKRead> locMap = new TreeMap<>();
        private Iterator<GATKRead> outputIterator = Collections.emptyIterator();
    }

    @Override
    protected void runTool( final JavaSparkContext ctx )
    {
        final JavaRDD<GATKRead> clusteredSplitReads = getReads().mapPartitions(
                readItr -> { return new SplitReadClusterer(readItr,getHeaderForReads().getSequenceDictionary()); } );

        try
        {
            ReadsSparkSink.writeReads(ctx, output, clusteredSplitReads, getHeaderForReads(), ReadsWriteFormat.SINGLE);
        }
        catch (IOException e)
        {
            throw new GATKException("unable to write bam: " + e);
        }
    }
}

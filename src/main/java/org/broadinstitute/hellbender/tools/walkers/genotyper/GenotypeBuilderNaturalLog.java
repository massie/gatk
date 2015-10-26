package org.broadinstitute.hellbender.tools.walkers.genotyper;

/**
 * This class is a hackish copy-and-paste of htsjdk.variant.variantcontext.GenotypeBuilder
 *
 * Unlike the original class, it assumes natural log likelihoods.  All changes are noted in the code.
 * Created by davidben on 10/25/15.
 */
import htsjdk.tribble.util.ParsingUtils;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.utils.QualityUtils;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A builder class for genotypes
 *
 * Provides convenience setter methods for all of the Genotype field
 * values.  Setter methods can be used in any order, allowing you to
 * pass through states that wouldn't be allowed in the highly regulated
 * immutable Genotype class.
 *
 * All fields default to meaningful MISSING values.
 *
 * Call make() to actually create the corresponding Genotype object from
 * this builder.  Can be called multiple times to create independent copies,
 * or with intervening sets to conveniently make similar Genotypes with
 * slight modifications.
 *
 * @author Mark DePristo
 * @since 06/12
 */
public final class GenotypeBuilderNaturalLog {
    private static final List<Allele> HAPLOID_NO_CALL = Arrays.asList(Allele.NO_CALL);
    private static final List<Allele> DIPLOID_NO_CALL = Arrays.asList(Allele.NO_CALL, Allele.NO_CALL);

    private String sampleName = null;
    private List<Allele> alleles = Collections.emptyList();

    private boolean isPhased = false;
    private int GQ = -1;
    private int DP = -1;
    private int[] AD = null;
    private int[] PL = null;
    private Map<String, Object> extendedAttributes = null;
    private String filters = null;
    private int initialAttributeMapSize = 5;

    private final static Map<String, Object> NO_ATTRIBUTES =
            Collections.unmodifiableMap(new HashMap<String, Object>(0));

    // -----------------------------------------------------------------
    //
    // Factory methods
    //
    // -----------------------------------------------------------------

    public static Genotype create(final String sampleName, final List<Allele> alleles) {
        return new GenotypeBuilderNaturalLog(sampleName, alleles).make();
    }

    public static Genotype create(final String sampleName,
                                  final List<Allele> alleles,
                                  final Map<String, Object> attributes) {
        return new GenotypeBuilderNaturalLog(sampleName, alleles).attributes(attributes).make();
    }

    protected static Genotype create(final String sampleName,
                                     final List<Allele> alleles,
                                     final double[] gls) {
        return new GenotypeBuilderNaturalLog(sampleName, alleles).PL(gls).make();
    }

    /**
     * Create a new Genotype object for a sample that's missing from the VC (i.e., in
     * the output header).  Defaults to a diploid no call genotype ./.
     *
     * @param sampleName the name of this sample
     * @return an initialized Genotype with sampleName that's a diploid ./. no call genotype
     */
    public static Genotype createMissing(final String sampleName, final int ploidy) {
        final GenotypeBuilderNaturalLog builder = new GenotypeBuilderNaturalLog(sampleName);
        switch ( ploidy ) {
            case 1:  builder.alleles(HAPLOID_NO_CALL); break;
            case 2:  builder.alleles(DIPLOID_NO_CALL); break;
            default: builder.alleles(Collections.nCopies(ploidy, Allele.NO_CALL)); break;
        }
        return builder.make();
    }

    /**
     * Create a empty builder.  Both a sampleName and alleles must be provided
     * before trying to make a Genotype from this builder.
     */
    public GenotypeBuilderNaturalLog() {}

    /**
     * Create a builder using sampleName.  Alleles must be provided
     * before trying to make a Genotype from this builder.
     * @param sampleName
     */
    public GenotypeBuilderNaturalLog(final String sampleName) {
        name(sampleName);
    }

    /**
     * Make a builder using sampleName and alleles for starting values
     * @param sampleName
     * @param alleles
     */
    public GenotypeBuilderNaturalLog(final String sampleName, final List<Allele> alleles) {
        name(sampleName);
        alleles(alleles);
    }

    /**
     * Create a new builder starting with the values in Genotype g
     * @param g
     */
    public GenotypeBuilderNaturalLog(final Genotype g) {
        copy(g);
    }

    /**
     * Copy all of the values for this builder from Genotype g
     * @param g
     * @return
     */
    public GenotypeBuilderNaturalLog copy(final Genotype g) {
        name(g.getSampleName());
        alleles(g.getAlleles());
        phased(g.isPhased());
        GQ(g.getGQ());
        DP(g.getDP());
        AD(g.getAD());
        PL(g.getPL());
        filter(g.getFilters());
        attributes(g.getExtendedAttributes());
        return this;
    }

    /**
     * Reset all of the builder attributes to their defaults.  After this
     * function you must provide sampleName and alleles before trying to
     * make more Genotypes.
     */
    public final void reset(final boolean keepSampleName) {
        if ( ! keepSampleName ) sampleName = null;
        alleles = Collections.emptyList();
        isPhased = false;
        GQ = -1;
        DP = -1;
        AD = null;
        PL = null;
        filters = null;
        extendedAttributes = null;
    }

    /**
     * Create a new Genotype object using the values set in this builder.
     *
     * After creation the values in this builder can be modified and more Genotypes
     * created, althrough the contents of array values like PL should never be modified
     * inline as they are not copied for efficiency reasons.
     *
     * @return a newly minted Genotype object with values provided from this builder
     */
    public Genotype make() {
        final Map<String, Object> ea = extendedAttributes == null ? NO_ATTRIBUTES : extendedAttributes;
        return new FastGenotypeNaturalLog(sampleName, alleles, isPhased, GQ, DP, AD, PL, filters, ea);
    }

    /**
     * Set this genotype's name
     * @param sampleName
     * @return
     */
    public GenotypeBuilderNaturalLog name(final String sampleName) {
        this.sampleName = sampleName;
        return this;
    }

    /**
     * Set this genotype's alleles
     * @param alleles
     * @return
     */
    public GenotypeBuilderNaturalLog alleles(final List<Allele> alleles) {
        if ( alleles == null )
            this.alleles = Collections.emptyList();
        else
            this.alleles = alleles;
        return this;
    }

    /**
     * Is this genotype phased?
     * @param phased
     * @return
     */
    public GenotypeBuilderNaturalLog phased(final boolean phased) {
        isPhased = phased;
        return this;
    }

    public GenotypeBuilderNaturalLog GQ(final int GQ) {
        this.GQ = GQ;
        return this;
    }

    /**  Set the GQ with a logPError value
     *
     * @param pLogError
     * @return
     */
    public GenotypeBuilderNaturalLog logPError(final double pLogError) {
        if ( pLogError == CommonInfo.NO_LOG10_PERROR )
            return noGQ();
        else
            return GQ((int)Math.round(QualityUtils.logProbToPhred(pLogError))); //CHANGED
    }

    /**
     * This genotype has no GQ value
     * @return
     */
    public GenotypeBuilderNaturalLog noGQ() { GQ = -1; return this; }

    /**
     * This genotype has no AD value
     * @return
     */
    public GenotypeBuilderNaturalLog noAD() { AD = null; return this; }

    /**
     * This genotype has no DP value
     * @return
     */
    public GenotypeBuilderNaturalLog noDP() { DP = -1; return this; }

    /**
     * This genotype has no PL value
     * @return
     */
    public GenotypeBuilderNaturalLog noPL() { PL = null; return this; }

    /**
     * This genotype has this DP value
     * @return
     */
    public GenotypeBuilderNaturalLog DP(final int DP) {
        this.DP = DP;
        return this;
    }

    /**
     * This genotype has this AD value
     * @return
     */
    public GenotypeBuilderNaturalLog AD(final int[] AD) {
        this.AD = AD;
        return this;
    }

    /**
     * This genotype has this PL value, as int[].  FAST
     * @return
     */
    public GenotypeBuilderNaturalLog PL(final int[] PL) {
        this.PL = PL;
        return this;
    }

    /**
     * This genotype has this PL value, converted from double[]. SLOW
     * @return
     */
    public GenotypeBuilderNaturalLog PL(final double[] GLs) {
        this.PL = GenotypeLikelihoodsNaturalLog.fromLogLikelihoods(GLs).getAsPLs();
        return this;
    }

    /**
     * This genotype has these attributes.
     *
     * Cannot contain inline attributes (DP, AD, GQ, PL)
     * @return
     */
    public GenotypeBuilderNaturalLog attributes(final Map<String, Object> attributes) {
        for ( Map.Entry<String, Object> pair : attributes.entrySet() )
            attribute(pair.getKey(), pair.getValue());
        return this;
    }

    /**
     * Tells this builder to remove all extended attributes
     *
     * @return
     */
    public GenotypeBuilderNaturalLog noAttributes() {
        this.extendedAttributes = null;
        return this;
    }

    /**
     * This genotype has this attribute key / value pair.
     *
     * Cannot contain inline attributes (DP, AD, GQ, PL)
     * @return
     */
    public GenotypeBuilderNaturalLog attribute(final String key, final Object value) {
        if ( extendedAttributes == null )
            extendedAttributes = new HashMap<String, Object>(initialAttributeMapSize);
        extendedAttributes.put(key, value);
        return this;
    }

    /**
     * Tells this builder to make a Genotype object that has had filters applied,
     * which may be empty (passes) or have some value indicating the reasons
     * why it's been filtered.
     *
     * @param filters non-null list of filters.  empty list =&gt; PASS
     * @return this builder
     */
    public GenotypeBuilderNaturalLog filters(final List<String> filters) {
        if ( filters.isEmpty() )
            return filter(null);
        else if ( filters.size() == 1 )
            return filter(filters.get(0));
        else
            return filter(ParsingUtils.join(";", ParsingUtils.sortList(filters)));
    }

    /**
     * varargs version of #filters
     * @param filters
     * @return
     */
    public GenotypeBuilderNaturalLog filters(final String ... filters) {
        return filters(Arrays.asList(filters));
    }

    /**
     * Most efficient version of setting filters -- just set the filters string to filters
     *
     * @param filter if filters == null or filters.equals("PASS") =&gt; genotype is PASS
     * @return
     */
    public GenotypeBuilderNaturalLog filter(final String filter) {
        this.filters = VCFConstants.PASSES_FILTERS_v4.equals(filter) ? null : filter;
        return this;
    }

    /**
     * This genotype is unfiltered
     *
     * @return
     */
    public GenotypeBuilderNaturalLog unfiltered() {
        return filter(null);
    }

    /**
     * Tell's this builder that we have at most these number of attributes
     * @return
     */
    public GenotypeBuilderNaturalLog maxAttributes(final int i) {
        initialAttributeMapSize = i;
        return this;
    }
}
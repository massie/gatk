package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.tribble.TribbleException;
import htsjdk.variant.utils.GeneralUtils;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.QualityUtils;

import java.util.Arrays;
import java.util.EnumMap;
import java.util.List;


/**
 * This class is a hackish copy-and-paste of htsjdk.variant.variantcontext.GenotypeLikelihoods
 *
 * Unlike the original class, it assumes natural log likelihoods.  All changes are noted in the code.
 * Created by davidben on 10/25/15.
 */
public class GenotypeLikelihoodsNaturalLog {


    private final static int NUM_LIKELIHOODS_CACHE_N_ALLELES = 5;
    private final static int NUM_LIKELIHOODS_CACHE_PLOIDY = 10;
    // caching numAlleles up to 5 and ploidy up to 10
    private final static int[][] numLikelihoodCache = new int[NUM_LIKELIHOODS_CACHE_N_ALLELES][NUM_LIKELIHOODS_CACHE_PLOIDY];

    public final static int MAX_PL = Integer.MAX_VALUE;

    //
    // There are two objects here because we are lazy in creating both representations
    // for this object: a vector of log10 Probs and the PL phred-scaled string.  Supports
    // having one set during initializating, and dynamic creation of the other, if needed
    //
    //CHANGE: this field renamed
    private double[] logLikelihoods = null;
    private String likelihoodsAsString_PLs = null;


    /**
     * initialize num likelihoods cache
     */
    static {
        // must be done before PLIndexToAlleleIndex
        for ( int numAlleles = 1; numAlleles < NUM_LIKELIHOODS_CACHE_N_ALLELES; numAlleles++ ) {
            for ( int ploidy = 1; ploidy < NUM_LIKELIHOODS_CACHE_PLOIDY; ploidy++ ) {
                numLikelihoodCache[numAlleles][ploidy] = calcNumLikelihoods(numAlleles, ploidy);
            }
        }
    }

    /**
     * The maximum number of alleles that we can represent as genotype likelihoods
     */
    public final static int MAX_ALT_ALLELES_THAT_CAN_BE_GENOTYPED = 50;

    /*
    * a cache of the PL index to the 2 alleles it represents over all possible numbers of alternate alleles
    */
    private final static GenotypeLikelihoodsNaturalLogAllelePair[] PLIndexToAlleleIndex = calculatePLcache(MAX_ALT_ALLELES_THAT_CAN_BE_GENOTYPED);

    public final static GenotypeLikelihoodsNaturalLog fromPLField(String PLs) {
        return new GenotypeLikelihoodsNaturalLog(PLs);
    }

    @Deprecated
    public final static GenotypeLikelihoodsNaturalLog fromGLField(String GLs) {
        return new GenotypeLikelihoodsNaturalLog(parseDeprecatedGLString(GLs));
    }

    //CHANGE: this method renamed
    public final static GenotypeLikelihoodsNaturalLog fromLogLikelihoods(double[] logLikelihoods) {
        return new GenotypeLikelihoodsNaturalLog(logLikelihoods);
    }

    //CHANGE: new method added for convenience
    public final static GenotypeLikelihoodsNaturalLog fromGenotype(final Genotype genotype) {
        return genotype.hasLikelihoods() ? fromPLs(genotype.getPL()) : null;
    }

    public final static GenotypeLikelihoodsNaturalLog fromPLs(final int[] pls) {
        return new GenotypeLikelihoodsNaturalLog(PLsToGLs(pls));
    }

    //
    // You must use the factory methods now
    //
    private GenotypeLikelihoodsNaturalLog(String asString) {
        likelihoodsAsString_PLs = asString;
    }

    private GenotypeLikelihoodsNaturalLog(double[] asVector) {
        logLikelihoods = asVector;
    }

    /**
     * Returns the genotypes likelihoods in negative log vector format.  pr{AA} = x, this
     * vector returns math.log10(x) for each of the genotypes.  Can return null if the
     * genotype likelihoods are "missing".
     *
     * @return
     */
    public double[] getAsVector() {
        // assumes one of the likelihoods vector or the string isn't null
        if ( logLikelihoods == null ) {
            // make sure we create the GL string if it doesn't already exist
            logLikelihoods = parsePLsIntoLikelihoods(likelihoodsAsString_PLs);
        }

        return logLikelihoods;
    }

    //CHANGED: new convenience methoid
    public static double[] likelihoodsFromGenotype(final Genotype genotype) {
        return genotype.hasLikelihoods() ? PLsToGLs(genotype.getPL()) : null;
    }


    public int[] getAsPLs() {
        final double[] GLs = getAsVector();
        return GLs == null ? null : GLsToPLs(GLs);
    }

    public String toString() {
        return getAsString();
    }

    public String getAsString() {
        if ( likelihoodsAsString_PLs == null ) {
            // todo -- should we accept null logLikelihoods and set PLs as MISSING?
            if ( logLikelihoods == null )
                throw new TribbleException("BUG: Attempted to get likelihoods as strings and neither the vector nor the string is set!");
            likelihoodsAsString_PLs = convertLikelihoodsToPLString(logLikelihoods);
        }

        return likelihoodsAsString_PLs;
    }

    @Override public boolean equals(Object aThat) {
        //check for self-comparison
        if ( this == aThat ) return true;

        if ( !(aThat instanceof GenotypeLikelihoodsNaturalLog) ) return false;
        GenotypeLikelihoodsNaturalLog that = (GenotypeLikelihoodsNaturalLog)aThat;

        // now a proper field-by-field evaluation can be made.
        // GLs are considered equal if the corresponding PLs are equal
        return Arrays.equals(getAsPLs(), that.getAsPLs());
    }

    //HACK -- THIS IS JUST OT REMOVE COMPILER WARNING!!!!
    @Override public int hashCode() {
        return logLikelihoods.hashCode();
    }

    //Return genotype likelihoods as an EnumMap with Genotypes as keys and likelihoods as values
    //Returns null in case of missing likelihoods
    //CHANGE: argument renamed
    public EnumMap<GenotypeType,Double> getAsMap(boolean normalizeFromLog){
        //Make sure that the loglikelihoods are set
        double[] likelihoods = normalizeFromLog ? GeneralUtils.normalizeFromLog10(getAsVector()) : getAsVector();
        if(likelihoods == null)
            return null;
        EnumMap<GenotypeType,Double> likelihoodsMap = new EnumMap<GenotypeType, Double>(GenotypeType.class);
        likelihoodsMap.put(GenotypeType.HOM_REF,likelihoods[GenotypeType.HOM_REF.ordinal()-1]);
        likelihoodsMap.put(GenotypeType.HET,likelihoods[GenotypeType.HET.ordinal()-1]);
        likelihoodsMap.put(GenotypeType.HOM_VAR, likelihoods[GenotypeType.HOM_VAR.ordinal() - 1]);
        return likelihoodsMap;
    }

    //Return the neg log Genotype Quality (GQ) for the given genotype
    //Returns Double.NEGATIVE_INFINITY in case of missing genotype

    /**
     * This is really dangerous and returns completely wrong results for genotypes from a multi-allelic context.
     * Use <code>getLogGQ(Genotype,VariantContext)</code>
     *  or <code>getLogGQ(Genotype,List&lt;Allele&gt;)</code> in place of it.
     *
     * If you <strong>know</strong> you're biallelic, use <code>getLogGQFromLikelihoods</code> directly.
     * @param genotype - actually a genotype type (no call, hom ref, het, hom var)
     * @return an unsafe quantity that could be negative. In the bi-allelic case, the GQ resulting from best minus next best (if the type is the best).
     */
    @Deprecated
    public double getLogGQ(GenotypeType genotype){
        return getLogGQFromLikelihoods(genotype.ordinal() - 1 /* NO_CALL IS FIRST */, getAsVector());
    }

    private double getLogGQ(List<Allele> genotypeAlleles, List<Allele> contextAlleles) {
        int allele1Index = contextAlleles.indexOf(genotypeAlleles.get(0));
        int allele2Index = contextAlleles.indexOf(genotypeAlleles.get(1));
        int plIndex = calculatePLindex(allele1Index,allele2Index);
        return getLogGQFromLikelihoods(plIndex, getAsVector());
    }

    public double getLogGQ(Genotype genotype, List<Allele> vcAlleles) {
        return getLogGQ(genotype.getAlleles(), vcAlleles);
    }

    public double getLogGQ(Genotype genotype, VariantContext context) {
        return getLogGQ(genotype, context.getAlleles());
    }

    public static double getLogGQFromLikelihoods(int iOfChoosenGenotype, double[] likelihoods){
        if(likelihoods == null)
            return Double.NEGATIVE_INFINITY;

        double qual = Double.NEGATIVE_INFINITY;
        for (int i=0; i < likelihoods.length; i++) {
            if (i==iOfChoosenGenotype)
                continue;
            if (likelihoods[i] >= qual)
                qual = likelihoods[i];
        }

        // qual contains now max(likelihoods[k]) for all k != bestGTguess
        qual = likelihoods[iOfChoosenGenotype] - qual;

        if (qual < 0) {
            // QUAL can be negative if the chosen genotype is not the most likely one individually.
            // In this case, we compute the actual genotype probability and QUAL is the likelihood of it not being the chosen one

            //CHANGE: this used to be normalization from log10 likelihoods
            //double[] normalized = GeneralUtils.normalizeFromLog10(likelihoods);   //OLD CODE
            double[] normalized = MathUtils.normalizeFromLog(likelihoods);          //NEW CODE
            double chosenGenotype = normalized[iOfChoosenGenotype];
            return Math.log(1.0 - chosenGenotype);                                //CHANGED TO NATURAL LOG
        } else {
            // invert the size, as this is the probability of making an error
            return -1 * qual;
        }
    }

    private final static double[] parsePLsIntoLikelihoods(String likelihoodsAsString_PLs) {
        if ( !likelihoodsAsString_PLs.equals(VCFConstants.MISSING_VALUE_v4) ) {
            String[] strings = likelihoodsAsString_PLs.split(",");
            double[] likelihoodsAsVector = new double[strings.length];
            try {
                for ( int i = 0; i < strings.length; i++ ) {
                    //CHANGED -- conversion factor from phred to log10 is -1/10
                    //conversion from phred to natural log is -log(10)/10
                    //likelihoodsAsVector[i] = Integer.parseInt(strings[i]) / -10.0;    //OLD CODE
                    likelihoodsAsVector[i] = QualityUtils.phredToLogProb(Integer.parseInt(strings[i]));      //NEW CODE
                }
            } catch (NumberFormatException e) {
                throw new TribbleException("The GL/PL tag contains non-integer values: " + likelihoodsAsString_PLs);
            }
            return likelihoodsAsVector;
        } else
            return null;
    }

    /**
     * Back-compatibility function to read old style GL formatted genotype likelihoods in VCF format
     * @param GLString
     * @return
     */
    private final static double[] parseDeprecatedGLString(String GLString) {
        if ( !GLString.equals(VCFConstants.MISSING_VALUE_v4) ) {
            String[] strings = GLString.split(",");
            double[] likelihoodsAsVector = new double[strings.length];
            for ( int i = 0; i < strings.length; i++ ) {
                likelihoodsAsVector[i] = Double.parseDouble(strings[i]);
            }
            return likelihoodsAsVector;
        }

        return null;
    }

    private final static String convertLikelihoodsToPLString(final double[] GLs) {
        if ( GLs == null )
            return VCFConstants.MISSING_VALUE_v4;

        final StringBuilder s = new StringBuilder();
        boolean first = true;
        for ( final int pl : GLsToPLs(GLs) ) {
            if ( ! first )
                s.append(",");
            else
                first = false;

            s.append(pl);
        }

        return s.toString();
    }

    private final static int[] GLsToPLs(final double[] GLs) {
        final int[] pls = new int[GLs.length];
        final double adjust = maxPL(GLs);

        for ( int i = 0; i < GLs.length; i++ ) {
            //CHANGED
            //pls[i] = (int)Math.round(Math.min(-10 * (GLs[i] - adjust), MAX_PL));  //OLD CODE
            pls[i] = (int)Math.round(Math.min(QualityUtils.logProbToPhred(GLs[i] - adjust), MAX_PL));   //NEW CODE
        }

        return pls;
    }

    private final static double maxPL(final double[] GLs) {
        double adjust = Double.NEGATIVE_INFINITY;
        for ( double l : GLs ) adjust = Math.max(adjust, l);
        return adjust;
    }

    private final static double[] PLsToGLs(final int pls[]) {
        double[] likelihoodsAsVector = new double[pls.length];
        for ( int i = 0; i < pls.length; i++ ) {
            //CHANGED
            //likelihoodsAsVector[i] = pls[i] / -10.0;                      //OLD CODE
            likelihoodsAsVector[i] = QualityUtils.phredToLogProb(pls[i]);   //NEW CODE
        }
        return likelihoodsAsVector;
    }

    // -------------------------------------------------------------------------------------
    //
    // Static conversion utilities, going from GL/PL index to allele index and vice versa.
    //
    // -------------------------------------------------------------------------------------

    /*
    * Class representing the 2 alleles (or rather their indexes into VariantContext.getAllele()) corresponding to a specific PL index.
    * Note that the reference allele is always index=0.
    */
    public static class GenotypeLikelihoodsNaturalLogAllelePair {
        public final int alleleIndex1, alleleIndex2;

        public GenotypeLikelihoodsNaturalLogAllelePair(final int alleleIndex1, final int alleleIndex2) {
            this.alleleIndex1 = alleleIndex1;
            this.alleleIndex2 = alleleIndex2;
        }
    }

    private static GenotypeLikelihoodsNaturalLogAllelePair[] calculatePLcache(final int altAlleles) {
        final int numLikelihoods = numLikelihoods(1 + altAlleles, 2);
        final GenotypeLikelihoodsNaturalLogAllelePair[] cache = new GenotypeLikelihoodsNaturalLogAllelePair[numLikelihoods];

        // for all possible combinations of 2 alleles
        for ( int allele1 = 0; allele1 <= altAlleles; allele1++ ) {
            for ( int allele2 = allele1; allele2 <= altAlleles; allele2++ ) {
                cache[calculatePLindex(allele1, allele2)] = new GenotypeLikelihoodsNaturalLogAllelePair(allele1, allele2);
            }
        }

        // a bit of sanity checking
        for ( int i = 0; i < cache.length; i++ ) {
            if ( cache[i] == null )
                throw new IllegalStateException("BUG: cache entry " + i + " is unexpected null");
        }

        return cache;
    }

    // -------------------------------------------------------------------------------------
    //
    // num likelihoods given number of alleles and ploidy
    //
    // -------------------------------------------------------------------------------------

    /**
     * Actually does the computation in @see #numLikelihoods
     *
     * @param numAlleles
     * @param ploidy
     * @return
     */
    private static final int calcNumLikelihoods(final int numAlleles, final int ploidy) {
        if (numAlleles == 1)
            return 1;
        else if (ploidy == 1)
            return numAlleles;
        else {
            int acc =0;
            for (int k=0; k <= ploidy; k++ )
                acc += calcNumLikelihoods(numAlleles - 1, ploidy - k);
            return acc;
        }
    }

    /**
     * Compute how many likelihood elements are associated with the given number of alleles
     * Equivalent to asking in how many ways N non-negative integers can add up to P is S(N,P)
     * where P = ploidy (number of chromosomes) and N = total # of alleles.
     * Each chromosome can be in one single state (0,...,N-1) and there are P of them.
     * Naive solution would be to store N*P likelihoods, but this is not necessary because we can't distinguish chromosome states, but rather
     * only total number of alt allele counts in all chromosomes.
     *
     * For example, S(3,2) = 6: For alleles A,B,C, on a diploid organism we have six possible genotypes:
     * AA,AB,BB,AC,BC,CC.
     * Another way of expressing is with vector (#of A alleles, # of B alleles, # of C alleles)
     * which is then, for ordering above, (2,0,0), (1,1,0), (0,2,0), (1,1,0), (0,1,1), (0,0,2)
     * In general, for P=2 (regular biallelic), then S(N,2) = N*(N+1)/2
     *
     * Note this method caches the value for most common num Allele / ploidy combinations for efficiency
     *
     * Recursive implementation:
     *   S(N,P) = sum_{k=0}^P S(N-1,P-k)
     *  because if we have N integers, we can condition 1 integer to be = k, and then N-1 integers have to sum to P-K
     * With initial conditions
     *   S(N,1) = N  (only way to have N integers add up to 1 is all-zeros except one element with a one. There are N of these vectors)
     *   S(1,P) = 1 (only way to have 1 integer add to P is with that integer P itself).
     *
     *   @param  numAlleles      Number of alleles (including ref)
     *   @param  ploidy          Ploidy, or number of chromosomes in set
     *   @return    Number of likelihood elements we need to hold.
     */
    public static int numLikelihoods(final int numAlleles, final int ploidy) {
        if ( numAlleles < NUM_LIKELIHOODS_CACHE_N_ALLELES
                && ploidy < NUM_LIKELIHOODS_CACHE_PLOIDY )
            return numLikelihoodCache[numAlleles][ploidy];
        else {
            // have to calculate on the fly
            return calcNumLikelihoods(numAlleles, ploidy);
        }
    }

    // As per the VCF spec: "the ordering of genotypes for the likelihoods is given by: F(j/k) = (k*(k+1)/2)+j.
    // In other words, for biallelic sites the ordering is: AA,AB,BB; for triallelic sites the ordering is: AA,AB,BB,AC,BC,CC, etc."
    // Assumes that allele1Index < allele2Index
    public static int calculatePLindex(final int allele1Index, final int allele2Index) {
        return (allele2Index * (allele2Index+1) / 2) + allele1Index;
    }

    /**
     * get the allele index pair for the given PL
     *
     * @param PLindex   the PL index
     * @return the allele index pair
     */
    public static GenotypeLikelihoodsNaturalLogAllelePair getAllelePair(final int PLindex) {
        // make sure that we've cached enough data
        if ( PLindex >= PLIndexToAlleleIndex.length )
            throw new IllegalStateException("Internal limitation: cannot genotype more than " + MAX_ALT_ALLELES_THAT_CAN_BE_GENOTYPED + " alleles");

        return PLIndexToAlleleIndex[PLindex];
    }

    // An index conversion from the deprecated PL ordering to the new VCF-based ordering for up to 3 alternate alleles
    protected static final int[] PLindexConversion = new int[]{0, 1, 3, 6, 2, 4, 7, 5, 8, 9};

    /**
     * get the allele index pair for the given PL using the deprecated PL ordering:
     *    AA,AB,AC,AD,BB,BC,BD,CC,CD,DD instead of AA,AB,BB,AC,BC,CC,AD,BD,CD,DD.
     * Although it's painful to keep this conversion around, our DiploidSNPGenotypeLikelihoods class uses the deprecated
     *    ordering and I know with certainty that external users have built code on top of it; changing it now would
     *    cause a whole lot of heartache for our collaborators, so for now at least there's a standard conversion method.
     * This method assumes at most 3 alternate alleles.
     *
     * @param PLindex   the PL index
     * @return the allele index pair
     */
    @Deprecated
    public static GenotypeLikelihoodsNaturalLogAllelePair getAllelePairUsingDeprecatedOrdering(final int PLindex) {
        return getAllelePair(PLindexConversion[PLindex]);
    }

    /**
     * get the PL indexes (AA, AB, BB) for the given allele pair; assumes allele1Index &lt;= allele2Index.
     *
     * @param allele1Index    the index in VariantContext.getAllele() of the first allele
     * @param allele2Index    the index in VariantContext.getAllele() of the second allele
     * @return the PL indexes
     */
    public static int[] getPLIndecesOfAlleles(final int allele1Index, final int allele2Index) {

        final int[] indexes = new int[3];
        indexes[0] = calculatePLindex(allele1Index, allele1Index);
        indexes[1] = calculatePLindex(allele1Index, allele2Index);
        indexes[2] = calculatePLindex(allele2Index, allele2Index);
        return indexes;
    }
}

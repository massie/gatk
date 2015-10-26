package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;

import java.util.List;
import java.util.Map;

/**
 * This is a horrible hack, just for getting the tests to pass before refactoring.
 *
 * It's literally identical to the htsjdk class except for making the construcotr public instead of protected
 * This is necessary for giving GenotypeBuilderNaturalLog access.
 * PERPETRATED by davidben on 10/26/15.
 */

/**
 * This class encompasses all the basic information about a genotype.  It is immutable.
 *
 * A genotype has several key fields
 *
 * <ul>
 * <li> a sample name, must be a non-null string</li>
 * <li> an ordered list of alleles, intrepreted as the genotype of the sample,
 *    each allele for each chromosome given in order.  If alleles = [a*, t]
 *    then the sample is a/t, with a (the reference from the *) the first
 *    chromosome and t on the second chromosome</li>
 * <li> an <code>isPhased</code> marker indicating where the alleles are phased with respect to some global
 *    coordinate system.  See VCF4.1 spec for a detailed discussion</li>
 * <li> Inline, optimized <code>int</code>s and <code>int[]</code> values for:
 * <ul>
 *      <li> GQ: the phred-scaled genotype quality, or <code>-1</code> if it's missing</li>
 *      <li> DP: the count of reads at this locus for this sample, or <code>-1</code> if missing</li>
 *      <li> AD: an array of counts of reads at this locus, one for each Allele at the site,
 *             that is, for each allele in the surrounding <code>VariantContext</code>.  <code>null</code> if missing.</li>
 *      <li> PL: phred-scaled genotype likelihoods in standard VCF4.1 order for
 *             all combinations of the alleles in the surrounding <code>VariantContext</code>, given
 *             the ploidy of the sample (from the alleles vector).  <code>null</code> if missing.</li>
 * </ul>
 * </li>
 *
 * <li> A general map from String keys to -&gt; Object values for all other attributes in
 *    this genotype.  Note that this map should not contain duplicate values for the
 *    standard bindings for GQ, DP, AD, and PL.  Genotype filters can be put into
 *    this genotype, but it isn't respected by the GATK in analyses</li>
 *</ul>
 *
 * <p>The only way to build a <code>Genotype</code> object is with a <code>GenotypeBuilder</code>, which permits values
 * to be set in any order, which means that <code>GenotypeBuilder</code> may at some in the chain of
 * sets pass through invalid states that are not permitted in a fully formed immutable
 * <code>Genotype</code>.</p>
 *
 * <p>Note this is a simplified, refactored Genotype object based on the original
 * generic (and slow) implementation from the original VariantContext + Genotype
 * codebase.</p>
 *
 * @author Mark DePristo
 * @since 05/12
 */
@SuppressWarnings("serial") //HORRID HACK ALERT!!!
public final class FastGenotypeNaturalLog extends Genotype {
    private final List<Allele> alleles;
    private final boolean isPhased;
    private final int GQ;
    private final int DP;
    private final int[] AD;
    private final int[] PL;
    private final Map<String, Object> extendedAttributes;

    /**
     * The only way to make one of these, for use by GenotypeBuilder only
     *
     * @param sampleName
     * @param alleles
     * @param isPhased
     * @param GQ
     * @param DP
     * @param AD
     * @param PL
     * @param extendedAttributes
     */
    public FastGenotypeNaturalLog(final String sampleName,
                           final List<Allele> alleles,
                           final boolean isPhased,
                           final int GQ,
                           final int DP,
                           final int[] AD,
                           final int[] PL,
                           final String filters,
                           final Map<String, Object> extendedAttributes) {
        super(sampleName, filters);
        this.alleles = alleles;
        this.isPhased = isPhased;
        this.GQ = GQ;
        this.DP = DP;
        this.AD = AD;
        this.PL = PL;
        this.extendedAttributes = extendedAttributes;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // Implmenting the abstract methods
    //
    // ---------------------------------------------------------------------------------------------------------

    @Override public List<Allele> getAlleles() {
        return alleles;
    }

    @Override public Allele getAllele(int i) {
        return alleles.get(i);
    }

    @Override public boolean isPhased() {
        return isPhased;
    }

    @Override public int getDP() {
        return DP;
    }

    @Override public int[] getAD() {
        return AD;
    }

    @Override public int getGQ()  {
        return GQ;
    }

    @Override public int[] getPL() {
        return PL;
    }

    // ---------------------------------------------------------------------------------------------------------
    //
    // get routines for extended attributes
    //
    // ---------------------------------------------------------------------------------------------------------

    public Map<String, Object> getExtendedAttributes() {
        return extendedAttributes;
    }

    /**
     * Is values a valid AD or PL field
     * @param values
     * @return
     */
    private static boolean validADorPLField(final int[] values) {
        if ( values != null )
            for ( int v : values )
                if ( v < 0 )
                    return false;
        return true;
    }
}
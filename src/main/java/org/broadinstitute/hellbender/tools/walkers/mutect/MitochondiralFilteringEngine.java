package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import org.broadinstitute.hellbender.utils.GATKProtectedVariantContextUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;

import java.util.Optional;

public class MitochondiralFilteringEngine extends Mutect2FilteringEngine {
    public MitochondiralFilteringEngine(MitochondrialFiltersArgumentCollection MTFAC, String tumorSample, Optional<String> normalSample, String lodKeyToUse) {
        super(MTFAC, tumorSample, normalSample, lodKeyToUse);
    }

    private void applyChimericOriginalAlignmentFilter(final MitochondrialFiltersArgumentCollection MTFAC, final VariantContext vc, final FilterResult filterResult) {
        final Genotype tumorGenotype = vc.getGenotype(tumorSample);
        final double[] alleleFractions = GATKProtectedVariantContextUtils.getAttributeAsDoubleArray(tumorGenotype, VCFConstants.ALLELE_FREQUENCY_KEY,
                () -> new double[] {1.0}, 1.0);
        final int maxFractionIndex = MathUtils.maxElementIndex(alleleFractions);
        final int[] ADs = tumorGenotype.getAD();
        final int altCount = ADs[maxFractionIndex + 1];

        if (tumorGenotype.hasAnyAttribute(GATKVCFConstants.ORIGINAL_CONTIG_MISMATCH_KEY) & vc.isBiallelic()) {
            int nonMtOa = Integer.parseInt(tumorGenotype.getAnyAttribute(GATKVCFConstants.ORIGINAL_CONTIG_MISMATCH_KEY).toString());
            if ((double) nonMtOa / altCount > MTFAC.nonMtAltByAlt) {
                filterResult.addFilter(GATKVCFConstants.CHIMERIC_ORIGINAL_ALIGNMENT_FILTER_NAME);
            }
        }
    }

    private void applyLODFilter(final MitochondrialFiltersArgumentCollection MTFAC, final VariantContext vc, final FilterResult filterResult) {
        if(vc.isBiallelic()) {
            Double LOD = vc.getAttributeAsDouble(GATKVCFConstants.LOD_KEY, 1);
            Double depth = vc.getAttributeAsDouble(VCFConstants.DEPTH_KEY, 1);
            Double lodByDepth = LOD / depth;
            if (lodByDepth < MTFAC.lodByDepth) {
                filterResult.addFilter(GATKVCFConstants.LOW_AVG_ALT_QUALITY_FILTER_NAME);
            }
        }
    }

    public FilterResult calculateFilters(MitochondrialFiltersArgumentCollection MTFAC, VariantContext vc) {
        final FilterResult filterResult = new FilterResult();

        applyInsufficientEvidenceFilter(MTFAC, vc, filterResult);
        applyDuplicatedAltReadFilter(MTFAC, vc, filterResult);
        applyStrandArtifactFilter(MTFAC, vc, filterResult);
        applyBaseQualityFilter(MTFAC, vc, filterResult);
        applyMappingQualityFilter(MTFAC, vc, filterResult);

        applyChimericOriginalAlignmentFilter(MTFAC, vc, filterResult);
        applyLODFilter(MTFAC, vc, filterResult);
        return filterResult;
    }
}

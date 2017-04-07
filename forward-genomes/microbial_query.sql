use Pfam;

# using the PfamA id
SELECT
    pfamseq.pfamseq_acc,
    swap.mbgd_genemap.sp_name,
    swap.mbgd_genemap.gene,
    pfamseq.auto_architecture,
    complete_proteomes.species,
    complete_proteomes.grouping,
    CONVERT( pfamseq.sequence USING UTF8),
    swap.mbgd_geneseq.sequence
FROM
    pfamseq
        INNER JOIN
    proteome_pfamseq ON pfamseq.auto_pfamseq = proteome_pfamseq.auto_pfamseq
        INNER JOIN
    complete_proteomes ON proteome_pfamseq.auto_proteome = complete_proteomes.auto_proteome
        INNER JOIN
    GeneCentric.topologyPreds ON pfamseq.pfamseq_acc = GeneCentric.topologyPreds.pfamseq_acc
        AND GeneCentric.topologyPreds.numSegments > 0
        AND GeneCentric.topologyPreds.predMethodKey = 1
        INNER JOIN
    swap.mbgd2pfam ON pfamseq.pfamseq_acc = swap.mbgd2pfam.pfamseq_acc
        INNER JOIN
    swap.mbgd_genemap ON swap.mbgd2pfam.sp_name = swap.mbgd_genemap.sp_name
        AND swap.mbgd2pfam.gene = swap.mbgd_genemap.gene
        INNER JOIN
    swap.mbgd_geneseq ON swap.mbgd2pfam.sp_name = swap.mbgd_geneseq.sp_name
        AND swap.mbgd2pfam.gene = swap.mbgd_geneseq.gene
WHERE
    complete_proteomes.auto_proteome IN (9, 36,
        84,
        369,
        444,
        695,
        753,
        885,
        1047,
        1241,
        1255,
        1323,
        1627,
        1707,
        1746,
        1792,
        1955,
        2169,
        2193,
        2291,
        2996,
        3025,
        3064)
        AND (LENGTH(pfamseq.sequence) + 1) * 3 = (swap.mbgd_genemap.to1 - swap.mbgd_genemap.from1 + 1)
GROUP BY pfamseq.pfamseq_acc
ORDER BY complete_proteomes.auto_proteome;

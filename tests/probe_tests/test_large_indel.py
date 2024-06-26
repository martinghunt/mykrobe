from base import assert_no_overlapping_kmers
from mykrobe.variants.schema.models import Variant


def test_large_variant1(variant_sets_and_reference, pg2):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="AACGCCCGGTATCTGAGGATCTGTGTTCTCACCCAATACAAGTCGCATTCACT",
        start=1355983,
        alternate_bases=["ACCGCCCGGTATCTGAGGATTGGTTTTCCACCCAAATACAAGTCGCATTCGCG"],
    )
    panel = pg2.create(v)
    assert_no_overlapping_kmers(panel)
    assert "TCGTCAACGCCCGGTATCTGAGGATCTGTGTTCTCACCCAATACAAGTCGCATTCACTGGA" in panel.refs
    assert panel.alts == [
        "TCGTCACCGCCCGGTATCTGAGGATTGGTTTTCCACCCAAATACAAGTCGCATTCGCGGGA"
    ]


def test_large_variant2(variant_sets_and_reference, pg2):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="AACGCCCGGTATCTGAGGATCTGTGTTCTCACCCAATACAAGTCGCATTCAC",
        start=1355983,
        alternate_bases=["ACCGCCCGGTATCTGAGGATTGGTTTTCCACCCAAATACAAGTCGCATTCGC"],
    )
    panel = pg2.create(v)
    assert_no_overlapping_kmers(panel)
    assert "TCGTCAACGCCCGGTATCTGAGGATCTGTGTTCTCACCCAATACAAGTCGCATTCACTGGA" in panel.refs
    assert panel.alts == [
        "TCGTCACCGCCCGGTATCTGAGGATTGGTTTTCCACCCAAATACAAGTCGCATTCGCTGGA"
    ]


def test_large_insertion(variant_sets_and_reference, pg2):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="C",
        start=2352065,
        alternate_bases=[
            "CCTCGCCTGGGCTGGCGAGCAGACGCAAAATCCCCCGCACGCCCGGCGTGTCGGGGGATTTTGCGTCTG"
        ],
    )
    panel = pg2.create(v)
    assert_no_overlapping_kmers(panel)
    assert "AGCTCGGCCAGCTCAGTCACGTCGCCGCCGCCTCGCCAGTTGACCGCGCCCGCTCGCGGCT" in panel.refs
    assert panel.alts == [
        "CCAGCTCAGTCACGTCGCCGCCGCCTCGCCTGGGCTGGCGAGCAGACGCAAAATCCCCCGCACGCCCGGCGTGTCGGGGGATTTTGCGTCTGCTCGCCAGTTGACCGCGCCCGCTCGCGGCT"
    ]


def test_large_var1(variant_sets_and_reference, pg2):
    variant_sets, reference = variant_sets_and_reference
    v = Variant.create(
        variant_sets=variant_sets,
        reference=reference,
        reference_bases="CGCGGGAGTAGAACGATCGCCAAGTGGTCGGTCTTGGCTGCCCACTTCATCCCCGGCGCCACCGGCAGGTCTCGCGGTCATCTCGACCAACGGAGGGCCGTCGGTGGTTCGTATCCGGCCAAGAACGGCGAGAACGGTTTGTGCCTCTATGCCAGGGTGAATGTCTCATCTCCCAGGCGGACGGTGATATCCAGTTCTCCGCCAAGAGCGGACACGTATTTGCGCAGTGTGTTGACCTGTGCGGAGCCGATGTCGCCGTTCTCGATGCTGGATACCCGGCTCTGCCGGATGTGCGCCAGCGCAGCCACCTGGACCTGGGTGAGTGACTGAGCCGCGCGCAGCTCCCGGAGCCGGAATGCCCGCACTTCATCGCGCATTCGTGCCTTGTGCCGGTCCACCGCCTCCCGGTTAACGGGACGTACGGCGTCCATGTCCCGTAGTGTCATCGCCATCGTGCCACTTACCCTTTCTTGCGCTTGCGCCTCTTTGGCTTCGTGTCCTCGAACTGTGCGAGATGTTCGGCAAACATCTCATCGGCCGCTTTGATCTTCTCGTCGTACCACTGGGTCCACCGCCCGGCCTTGTTACCGGCGGCCAGCATGATCGCCTGCCGCGCCGGGTCGAAGGCGAACAGAATGCGGACCTCGGACCGCCCTTGTGATCCTGGACGCAGCTCCTTCATGTTCTTGTGGCGCGACCCACGCACCGTGTCCACCAGAGGACAGCCAAGTGCGGGGCCCTCTTCCTCGAGAACCTCGATAGCTGCGAACACCAATTCGTAGGTCTCTCGGTCCAAGCCGTTGAGCCAGGCGGAGATGCGCTCCACATCCGCCGTCCACCCCACAGAGTCGCAGAGTAGCGCGATACGCGATATCACACAAGGGTGATATTCCTCCGGGTAAGAGCAGCGGGCGACGGGGCTACCGTCGAGGAAATGCCGGCAGGCGAGGACGGACTCTGCGCACCCGGGCCGTTGAAACAGTAGCCTGTGCCAGGCCGAGAATTCATCCCCACGTATGAGGCAGTACAGTGCGCCGCCGTGCGCGTTCTCCCATGGAACGTTCACGGGCTCCCGTGGATGACAGGCGTTTCATGAACGCCAGCGCCGCCGCAACCCGACCGAAAGCGGTTGACCCCAAGGAGAGCTGGAAGTCGAGGCCACCACCTTCGCCGCGGAGTTGCTCATGCCCGAGAGCGAGACTCGTCCCGAAATACGCCGGCTCGATTTCGGCAAGTTGCTCGAACTGAAGCGGGAATGGGCGTCGACCCGCTCGACCAGCCCCAGCCGGGTGACCAGCCCCAGCCGGGTGACCAGCCGATGCACCGCGGCGATCCCACCGAAGCCGGTGGCATCGATGTTGGCGCCGACCTCGTAGCGCACCGCGCCCGAACCCAGCATCGGCCTGGGCTGCGCCGCCCAGCGTCCAGCCCGCGCGTGCCGCGCCGCCACCCTGCGCCCTCGGCGTGTGATGTTTCGCCGACTCTGTTCATGGGTTATCTTCTTCACCACAAAGGCCTTTCCTGCTGGGCTGTGTTGAGGTCGCAAACCCAGCCAGGGTAAGGCCTTTGGCCTCTCCTACCCGGCCGACACGCTTACTGAAGGCCTAGTCTAGGCAGGCCATTCAATCTGCGGAATCGAAAAATTCGGTTCCAGCCTGCTCGTTTCCTTTCCGACAGCGATCTGACGTTGCGTAACGTCATTTGTACGGACTCTTTTAGCGGCATTGATTTCAGATGCCAACGCCGTCTGTGCTGTAGCGCCGATTGGCCGAAACTGTAAATTTGTATGATTATTTAAATCTTTGACGAACACGCGCCACAAACGTACTATCTCTTTGGCAAAGTCCACCGGCATCTCATTCAACGGTTTTGTTTGCGCGTGGTCGTCATATGTTGGTAACTGTGTAACCGGCCGCCTATCTTGCGCGTGCATCATATGACTATGAATCGGCCTTCTCCAGTGAAATTGATACAAGATCGATCCGATAAGCGGTACCTTGTACACAGTGCAATTGTAGTAATTCGCGTTTTGTCCTACGCTTGTATTCTGCGTGAAGAATTCA",
        start=2266659,
        alternate_bases=[
            "CACGCGAGTTGTAGATGATCGTTGAGTGGTCTTGCTTGGACTTCCATTTCATCTTTTCGACGCGCCAGGTCTCGCGGTCCTCCGGATCTGCGCCCGGTTTGAGTTGCACATCAAGGGGATACGGCTTGACCGACTCGTAGCCGACATGTAAGTCGGCTAGTTTCCGGCCGGCGCTGGCGAGCTGGTCGAAGCGTTCGCGGGTCTCCGGTGTTGGGATGTGCGGGAGCATCTTCTTGAGGTCAGCGGCGTATTTTGTGCGGTAGGCGGGGTCATGCAGCAGGCCGTAGACGTAGTAGAAGATGTCGTCTTTGGTGACTTGGTCGCCGATCGTGTCGCGGTAGAGCTTGAGGATGACGCCGGTGATGTTGTCGACGCGGCGGTAGCCGTGGTCGTCTACTTCGGCGTTGGTGGTGGACTCGAAATCGAGTTCGCCGTCACGTGGTTCGGTCTTCTCGTAGGTCCAGCGCGGGAAGAATTGACCGTTGCTTGAGCCCCAGAATGCGAGATCGGGGATAGCGTTTAGCATCAGACACGAGAAGGGCTTGTCTGAGCCCATGCCAACCACGTAGTAACCGACATTCCCGTGCTCCGGCGTCGGAAACATCGACGGAAGCTGGTAGGTACAGTTGTTGAGCTGCTGGTTGGGGTCGAGGTAGGCGTGCTCTTTCGTAAATGGTCGGTACGTGCCGAGCCGCATTCCCGCGGGAGCGAATTCGATGCGAATGCCTTGTGCCACTTGCCGCTTGTTGATGCGGTCCCAGCTGAACTTGGCCGAGTCCACGGTAATGAGGGCGTCAACCGGCGGGGTCTTGGCGTCCCTTCCGCGGATCTCGTTGATCCGGTCGACCTCCGAGTTGTAGAAGTCGATCGTGCGTCCGATGTTGGCCTCGAGCGCACCACGTGAAAAGTTGTAACACCACGCATCCCGGCTGGTCTTCAAGCCCGCGGAATAGTTCGCGAAGACACGTGTCACGTCAAGAGCAGCCTTCTTGTCGCCGATAACCGGCCACGCGCTGAACGCGTCGTCGCGTTGGTTGACCCAGTCACCGTGCAAGTTGGGTGTGACTGTCTGCCATTCCACCGTGTCGAGGTAGCCGTCGCCGACGATCCGCAACTTCTCCTCGCGACTCAGGTAATCGCCGATGTCGCGGTAAAGGACATCGCATGGCCCGCTGTGCTTCGGATCCTTGATGCCAAGGAAGATCGCCACCGTGTTGCGACTCCCCCCGCCAAAGACCTTGCCGCCTTCCTGGCGTGAGAGTTCCCCAGCTGTGCGCTGGTTCCCCCGCAGGTTGTACACATATACCGCCGCGTAGTCGTCGGCGAGCGACAACCGCATGCCGTCTGCCGTGTTGCCGTCTATGTACCCACCATTGGAGACGAATCCGACAACACCGTTGTCACCAATGCGGTCGGTCGCCCACCGGAACGCGCGAATATACGAGTCGTACAGGCTGTTCTTCAGCTGCGCCGTCGACCGCTTCGCGTACGTCTGCTCAATCCGCCCGTCCAACGTCGGATACTTCACGTTGGCGTTCAGGTCGTTCGCGCTGCTCTGCCCCACCGAGTACGGCGGATTCCCGATGATCACGCTGATCGGCGTCGCCAGCTGTCGCAAGATCCGAGCGTTGTTGTACGGGAACATGATCGCGTCCATCGAGTCCCCGGCTTCGGAAATCTGGAACGTGTCGGCCAGCGCCATCCCGGGGAACGGCTCATAGGCGTCGGCGTCGGCGGTCTTGCCCGCCAAAGCATGGTAGGTCGACTCGATGTTCACCGCGGCGATGTAGTACGCCAGCAGCATGATCTCGTTGGCGTGCAGCTCTTGCGAGTACTTTCGGGTGAGGTCGGCGGCCGTGATCAGGTCGGACTGCAGCAGCCGGGTAATGAATGTGCCCGTCCCGGCGAAGCCGTCCAGAATATGCACGCCCTCGTCGGTCAGCCCGCGCCCGAAATGCTTGCGCGACACGAAATCAGCCGCCCGCACAATGAAGTCCACGACCTCGACCGGCGTGTACACGATCCCCAGCGCCTCGGCCTGCTTCTTGAAGCCGATGCGGAAGAACTTCTCGTACAGCTCGGCGATCACCTGCTGCTTGCCCTCGGCGCTGGTGACCTCGCCGGCGCGCCGTCGCACCGATTCGTAAAAGCCTTCCAACCGAGCGGTTTCGGCCTCCAGGCCGGCACCCCCGACGGTGTCGACCATCTTCTGCATGGCCCGCGACACCGGGTTGTGCGACGCGAAGTCATGCCCGGCGAACAGCGCGTCGAACACCGGCTTGGTGATCAGGTGCTGCGAGAGCATGCTGATCGCGTCATCGGGGGTGATCGAGTCATTGAGGTTATCGCGCAGCCCGGCCAGGAACTGCTCGAACGCCGCCGCCGCCGTAGCGTCGGCGCCGCCGAGCAGGGCGTGGATACGGGTGGTCAGCGTCGCGGCGATGTCGGCGACATCGGCGGCCCACTGCTCCCAATAGGTCCGGGTGCCAACCTTGTCGACGATGCGCGCGTAGATCGCTTCCTGCCACTGCGACAACGAGAACATCGCCAACTGCTCCGCGACGGCGGGTCCCGCCTCGTCGGAGGTCGGCCCGATGTGACCGCCCAACAGCTTGTCGCTGCCTTCACCGGTCTTCGTCGGCTTCACGTTCAGCGCAATGCTGTTCACCATCGCGTCGAAGCGCTCGTCGTGCGACCGCAACGCGTTGAGGACCTGCCACACCACCTTGAACCGTTTGTTGTCGGCCAACGCGGCAGACGGCTCGACACCCTCGGGCACCGCCACCGGCAAGATGACGTACCCGTAGTCCTTGCCGGGCGACTTGCGCATCACCCGACCGACCGACTGCACCACGTCGACGATGGAATTGCGCGGATTCAGGAACAGCACCGCGTCCAGCGCGGGCACGTCGACCCCTTCGGAGAGGCAGCGGGCGTTGGACAGGATGCGGCATTCATCCTCGGCGACCACGCCTTTGAGCCAGGCCAGCTGTTCGTTGCGGACCAGCGCGTTGAACGTCCCGTCCACGTGGCGCACCG"
        ],
    )
    panel = pg2.create(v)
    assert_no_overlapping_kmers(panel)
    assert (
        "TGGTGACGCGGGAGTAGAACGATCGCCAAGTGGTCGGTCTTGGCTGCCCACTTCATCCCCGGCGCCACCGGCAGGTCTCGCGGTCATCTCGACCAACGGAGGGCCGTCGGTGGTTCGTATCCGGCCAAGAACGGCGAGAACGGTTTGTGCCTCTATGCCAGGGTGAATGTCTCATCTCCCAGGCGGACGGTGATATCCAGTTCTCCGCCAAGAGCGGACACGTATTTGCGCAGTGTGTTGACCTGTGCGGAGCCGATGTCGCCGTTCTCGATGCTGGATACCCGGCTCTGCCGGATGTGCGCCAGCGCAGCCACCTGGACCTGGGTGAGTGACTGAGCCGCGCGCAGCTCCCGGAGCCGGAATGCCCGCACTTCATCGCGCATTCGTGCCTTGTGCCGGTCCACCGCCTCCCGGTTAACGGGACGTACGGCGTCCATGTCCCGTAGTGTCATCGCCATCGTGCCACTTACCCTTTCTTGCGCTTGCGCCTCTTTGGCTTCGTGTCCTCGAACTGTGCGAGATGTTCGGCAAACATCTCATCGGCCGCTTTGATCTTCTCGTCGTACCACTGGGTCCACCGCCCGGCCTTGTTACCGGCGGCCAGCATGATCGCCTGCCGCGCCGGGTCGAAGGCGAACAGAATGCGGACCTCGGACCGCCCTTGTGATCCTGGACGCAGCTCCTTCATGTTCTTGTGGCGCGACCCACGCACCGTGTCCACCAGAGGACAGCCAAGTGCGGGGCCCTCTTCCTCGAGAACCTCGATAGCTGCGAACACCAATTCGTAGGTCTCTCGGTCCAAGCCGTTGAGCCAGGCGGAGATGCGCTCCACATCCGCCGTCCACCCCACAGAGTCGCAGAGTAGCGCGATACGCGATATCACACAAGGGTGATATTCCTCCGGGTAAGAGCAGCGGGCGACGGGGCTACCGTCGAGGAAATGCCGGCAGGCGAGGACGGACTCTGCGCACCCGGGCCGTTGAAACAGTAGCCTGTGCCAGGCCGAGAATTCATCCCCACGTATGAGGCAGTACAGTGCGCCGCCGTGCGCGTTCTCCCATGGAACGTTCACGGGCTCCCGTGGATGACAGGCGTTTCATGAACGCCAGCGCCGCCGCAACCCGACCGAAAGCGGTTGACCCCAAGGAGAGCTGGAAGTCGAGGCCACCACCTTCGCCGCGGAGTTGCTCATGCCCGAGAGCGAGACTCGTCCCGAAATACGCCGGCTCGATTTCGGCAAGTTGCTCGAACTGAAGCGGGAATGGGCGTCGACCCGCTCGACCAGCCCCAGCCGGGTGACCAGCCCCAGCCGGGTGACCAGCCGATGCACCGCGGCGATCCCACCGAAGCCGGTGGCATCGATGTTGGCGCCGACCTCGTAGCGCACCGCGCCCGAACCCAGCATCGGCCTGGGCTGCGCCGCCCAGCGTCCAGCCCGCGCGTGCCGCGCCGCCACCCTGCGCCCTCGGCGTGTGATGTTTCGCCGACTCTGTTCATGGGTTATCTTCTTCACCACAAAGGCCTTTCCTGCTGGGCTGTGTTGAGGTCGCAAACCCAGCCAGGGTAAGGCCTTTGGCCTCTCCTACCCGGCCGACACGCTTACTGAAGGCCTAGTCTAGGCAGGCCATTCAATCTGCGGAATCGAAAAATTCGGTTCCAGCCTGCTCGTTTCCTTTCCGACAGCGATCTGACGTTGCGTAACGTCATTTGTACGGACTCTTTTAGCGGCATTGATTTCAGATGCCAACGCCGTCTGTGCTGTAGCGCCGATTGGCCGAAACTGTAAATTTGTATGATTATTTAAATCTTTGACGAACACGCGCCACAAACGTACTATCTCTTTGGCAAAGTCCACCGGCATCTCATTCAACGGTTTTGTTTGCGCGTGGTCGTCATATGTTGGTAACTGTGTAACCGGCCGCCTATCTTGCGCGTGCATCATATGACTATGAATCGGCCTTCTCCAGTGAAATTGATACAAGATCGATCCGATAAGCGGTACCTTGTACACAGTGCAATTGTAGTAATTCGCGTTTTGTCCTACGCTTGTATTCTGCGTGAAGAATTCAAACA"
        in panel.refs
    )
    assert panel.alts == [
        "GGCCTCGTCGGGAATGCCGGCGATGGTGACACGCGAGTTGTAGATGATCGTTGAGTGGTCTTGCTTGGACTTCCATTTCATCTTTTCGACGCGCCAGGTCTCGCGGTCCTCCGGATCTGCGCCCGGTTTGAGTTGCACATCAAGGGGATACGGCTTGACCGACTCGTAGCCGACATGTAAGTCGGCTAGTTTCCGGCCGGCGCTGGCGAGCTGGTCGAAGCGTTCGCGGGTCTCCGGTGTTGGGATGTGCGGGAGCATCTTCTTGAGGTCAGCGGCGTATTTTGTGCGGTAGGCGGGGTCATGCAGCAGGCCGTAGACGTAGTAGAAGATGTCGTCTTTGGTGACTTGGTCGCCGATCGTGTCGCGGTAGAGCTTGAGGATGACGCCGGTGATGTTGTCGACGCGGCGGTAGCCGTGGTCGTCTACTTCGGCGTTGGTGGTGGACTCGAAATCGAGTTCGCCGTCACGTGGTTCGGTCTTCTCGTAGGTCCAGCGCGGGAAGAATTGACCGTTGCTTGAGCCCCAGAATGCGAGATCGGGGATAGCGTTTAGCATCAGACACGAGAAGGGCTTGTCTGAGCCCATGCCAACCACGTAGTAACCGACATTCCCGTGCTCCGGCGTCGGAAACATCGACGGAAGCTGGTAGGTACAGTTGTTGAGCTGCTGGTTGGGGTCGAGGTAGGCGTGCTCTTTCGTAAATGGTCGGTACGTGCCGAGCCGCATTCCCGCGGGAGCGAATTCGATGCGAATGCCTTGTGCCACTTGCCGCTTGTTGATGCGGTCCCAGCTGAACTTGGCCGAGTCCACGGTAATGAGGGCGTCAACCGGCGGGGTCTTGGCGTCCCTTCCGCGGATCTCGTTGATCCGGTCGACCTCCGAGTTGTAGAAGTCGATCGTGCGTCCGATGTTGGCCTCGAGCGCACCACGTGAAAAGTTGTAACACCACGCATCCCGGCTGGTCTTCAAGCCCGCGGAATAGTTCGCGAAGACACGTGTCACGTCAAGAGCAGCCTTCTTGTCGCCGATAACCGGCCACGCGCTGAACGCGTCGTCGCGTTGGTTGACCCAGTCACCGTGCAAGTTGGGTGTGACTGTCTGCCATTCCACCGTGTCGAGGTAGCCGTCGCCGACGATCCGCAACTTCTCCTCGCGACTCAGGTAATCGCCGATGTCGCGGTAAAGGACATCGCATGGCCCGCTGTGCTTCGGATCCTTGATGCCAAGGAAGATCGCCACCGTGTTGCGACTCCCCCCGCCAAAGACCTTGCCGCCTTCCTGGCGTGAGAGTTCCCCAGCTGTGCGCTGGTTCCCCCGCAGGTTGTACACATATACCGCCGCGTAGTCGTCGGCGAGCGACAACCGCATGCCGTCTGCCGTGTTGCCGTCTATGTACCCACCATTGGAGACGAATCCGACAACACCGTTGTCACCAATGCGGTCGGTCGCCCACCGGAACGCGCGAATATACGAGTCGTACAGGCTGTTCTTCAGCTGCGCCGTCGACCGCTTCGCGTACGTCTGCTCAATCCGCCCGTCCAACGTCGGATACTTCACGTTGGCGTTCAGGTCGTTCGCGCTGCTCTGCCCCACCGAGTACGGCGGATTCCCGATGATCACGCTGATCGGCGTCGCCAGCTGTCGCAAGATCCGAGCGTTGTTGTACGGGAACATGATCGCGTCCATCGAGTCCCCGGCTTCGGAAATCTGGAACGTGTCGGCCAGCGCCATCCCGGGGAACGGCTCATAGGCGTCGGCGTCGGCGGTCTTGCCCGCCAAAGCATGGTAGGTCGACTCGATGTTCACCGCGGCGATGTAGTACGCCAGCAGCATGATCTCGTTGGCGTGCAGCTCTTGCGAGTACTTTCGGGTGAGGTCGGCGGCCGTGATCAGGTCGGACTGCAGCAGCCGGGTAATGAATGTGCCCGTCCCGGCGAAGCCGTCCAGAATATGCACGCCCTCGTCGGTCAGCCCGCGCCCGAAATGCTTGCGCGACACGAAATCAGCCGCCCGCACAATGAAGTCCACGACCTCGACCGGCGTGTACACGATCCCCAGCGCCTCGGCCTGCTTCTTGAAGCCGATGCGGAAGAACTTCTCGTACAGCTCGGCGATCACCTGCTGCTTGCCCTCGGCGCTGGTGACCTCGCCGGCGCGCCGTCGCACCGATTCGTAAAAGCCTTCCAACCGAGCGGTTTCGGCCTCCAGGCCGGCACCCCCGACGGTGTCGACCATCTTCTGCATGGCCCGCGACACCGGGTTGTGCGACGCGAAGTCATGCCCGGCGAACAGCGCGTCGAACACCGGCTTGGTGATCAGGTGCTGCGAGAGCATGCTGATCGCGTCATCGGGGGTGATCGAGTCATTGAGGTTATCGCGCAGCCCGGCCAGGAACTGCTCGAACGCCGCCGCCGCCGTAGCGTCGGCGCCGCCGAGCAGGGCGTGGATACGGGTGGTCAGCGTCGCGGCGATGTCGGCGACATCGGCGGCCCACTGCTCCCAATAGGTCCGGGTGCCAACCTTGTCGACGATGCGCGCGTAGATCGCTTCCTGCCACTGCGACAACGAGAACATCGCCAACTGCTCCGCGACGGCGGGTCCCGCCTCGTCGGAGGTCGGCCCGATGTGACCGCCCAACAGCTTGTCGCTGCCTTCACCGGTCTTCGTCGGCTTCACGTTCAGCGCAATGCTGTTCACCATCGCGTCGAAGCGCTCGTCGTGCGACCGCAACGCGTTGAGGACCTGCCACACCACCTTGAACCGTTTGTTGTCGGCCAACGCGGCAGACGGCTCGACACCCTCGGGCACCGCCACCGGCAAGATGACGTACCCGTAGTCCTTGCCGGGCGACTTGCGCATCACCCGACCGACCGACTGCACCACGTCGACGATGGAATTGCGCGGATTCAGGAACAGCACCGCGTCCAGCGCGGGCACGTCGACCCCTTCGGAGAGGCAGCGGGCGTTGGACAGGATGCGGCATTCATCCTCGGCGACCACGCCTTTGAGCCAGGCCAGCTGTTCGTTGCGGACCAGCGCGTTGAACGTCCCGTCCACGTGGCGCACCGAACACGCCAGGCCCGGGCCGTCGTCAACCA"
    ]

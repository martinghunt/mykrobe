from unittest import TestCase
from mykrobe.variants.schema.models import Variant
from mykrobe.variants.schema.models import VariantCall
from mykrobe.typing import VariantTyper
from mykrobe.typing import ProbeCoverage
from mykrobe.typing import SequenceProbeCoverage
from mykrobe.typing import VariantProbeCoverage


class VariantTyperTest(TestCase):

    def setUp(self):
        self.vt = VariantTyper(expected_depths=[100])

    def teardown(self):
        pass

    def test_wt_vars(self):
        reference_coverage = ProbeCoverage(min_depth=100,
                                           percent_coverage=100,
                                           median_depth=100,
                                           k_count=100)
        alternate_coverages = [ProbeCoverage(min_depth=100,
                                             percent_coverage=3,
                                             median_depth=100,
                                             k_count=3)]
        v1 = VariantProbeCoverage(var_name="A123T",
                                  reference_coverages=[reference_coverage],
                                  alternate_coverages=alternate_coverages
                                  )
        call = self.vt.type([v1])
        self.assertEqual([0, 0], call['genotype'])
        self.assertEqual([100], call["info"]['expected_depths'])

    def test_alt_vars(self):
        reference_coverage = ProbeCoverage(min_depth=100,
                                           percent_coverage=3,
                                           median_depth=100,
                                           k_count=3)
        alternate_coverages = [ProbeCoverage(min_depth=100,
                                             percent_coverage=100,
                                             median_depth=100,
                                             k_count=100)]
        v1 = VariantProbeCoverage(var_name="A123T",
                                  reference_coverages=[reference_coverage],
                                  alternate_coverages=alternate_coverages
                                  )
        call = self.vt.type([v1])
        self.assertEqual([1, 1], call['genotype'])

    def test_mixed_vars(self):
        reference_coverage = ProbeCoverage(min_depth=100,
                                           percent_coverage=100,
                                           median_depth=50,
                                           k_count=50)
        alternate_coverages = [ProbeCoverage(min_depth=100,
                                             percent_coverage=100,
                                             median_depth=50,
                                             k_count=50)]
        v1 = VariantProbeCoverage(var_name="A123T",
                                  reference_coverages=[reference_coverage],
                                  alternate_coverages=alternate_coverages
                                  )
        call = self.vt.type(v1)
        self.assertEqual([0, 1], call['genotype'])

    def test_mixed_vars2(self):
        reference_coverage = ProbeCoverage(min_depth=11,
                                           percent_coverage=100,
                                           median_depth=42,
                                           k_count=42)
        alternate_coverages = [ProbeCoverage(min_depth=94,
                                             percent_coverage=100,
                                             median_depth=102,
                                             k_count=94)]
        v1 = VariantProbeCoverage(var_name="A123T",
                                  reference_coverages=[reference_coverage],
                                  alternate_coverages=alternate_coverages
                                  )
        call = self.vt.type(v1)
        self.assertEqual([0, 1], call['genotype'])


class VariantTyperWithContamination(TestCase):

    def setUp(self):
        self.vt_no_contaim = VariantTyper(
            expected_depths=[100],
            contamination_depths=[])
        # To do add contamination type
        # self.vt_contaim = VariantTyper(
        #     expected_depths=[80],
        #     contamination_depths=[20])

    def teardown(self):
        pass

    def test_simple_case(self):
        reference_coverage = ProbeCoverage(min_depth=100,
                                           percent_coverage=100,
                                           median_depth=80,
                                           k_count=80)
        alternate_coverages = [ProbeCoverage(min_depth=100,
                                             percent_coverage=100,
                                             median_depth=20,
                                             k_count=40)]
        v1 = VariantProbeCoverage(var_name="A123T",
                                  reference_coverages=[reference_coverage],
                                  alternate_coverages=alternate_coverages
                                  )

        call = self.vt_no_contaim.type(v1)
        self.assertEqual([0, 1], call['genotype'])

        # call = self.vt_contaim.type(v1)
        # assert call['genotype'] == [0, 0]


class TestVariantTyperWithMultipleAlternateCoverages(TestCase):

    def setUp(self):
        # to do, test should pass on kc model also
        self.vt_no_contaim = VariantTyper(
            expected_depths=[100],
            contamination_depths=[],
            model="median_depth")

    def teardown(self):
        pass

    def test_simple_case(self):
        reference_coverage = ProbeCoverage(min_depth=100,
                                           percent_coverage=70,
                                           median_depth=80,
                                           k_count=80)
        alt1 = ProbeCoverage(min_depth=100,
                             percent_coverage=70,
                             median_depth=20,
                             k_count=20)
        alt2 = ProbeCoverage(min_depth=100,
                             percent_coverage=100,
                             median_depth=80,
                             k_count=80)
        alternate_coverages = [alt1, alt2]
        v1 = VariantProbeCoverage(var_name="A123T",
                                  reference_coverages=[reference_coverage],
                                  alternate_coverages=alternate_coverages
                                  )
        self.assertEqual(alt2, v1._choose_best_alternate_coverage())

        call = self.vt_no_contaim.type(v1)
        self.assertEqual([1, 1], call['genotype'])


class TestVariantTyperWithMultipleProbeCoverages(TestCase):

    def setUp(self):
        self.vt_no_contaim = VariantTyper(
            expected_depths=[100],
            contamination_depths=[])

    def teardown(self):
        pass

    def test_simple_case(self):
        reference_coverage = ProbeCoverage(min_depth=100,
                                           percent_coverage=80,
                                           median_depth=80,
                                           k_count=80)
        alt1 = ProbeCoverage(min_depth=100,
                             percent_coverage=50,
                             median_depth=20,
                             k_count=20)
        alt2 = ProbeCoverage(min_depth=100,
                             percent_coverage=40,
                             median_depth=80,
                             k_count=30)
        alternate_coverages = [alt1, alt2]

        v1 = VariantProbeCoverage(var_name="A123T",
                                  reference_coverages=[reference_coverage],
                                  alternate_coverages=alternate_coverages
                                  )

        reference_coverage = ProbeCoverage(min_depth=100,
                                           percent_coverage=80,
                                           median_depth=80,
                                           k_count=20)
        alt1 = ProbeCoverage(min_depth=100,
                             percent_coverage=50,
                             median_depth=20,
                             k_count=20)
        alt2 = ProbeCoverage(min_depth=100,
                             percent_coverage=100,
                             median_depth=80,
                             k_count=100)

        alternate_coverages = [alt1, alt2]

        v2 = VariantProbeCoverage(var_name="A123T",
                                  reference_coverages=[reference_coverage],
                                  alternate_coverages=alternate_coverages
                                  )

        call = self.vt_no_contaim.type([v1, v2])
        self.assertEqual([1, 1], call['genotype'])


class TestVariantTyperWithLowMinimum(TestCase):

    def setUp(self):
        self.vt_no_contaim = VariantTyper(
            expected_depths=[100],
            contamination_depths=[])
        self.vt2_no_contaim = VariantTyper(
            expected_depths=[1],
            contamination_depths=[])

    def teardown(self):
        pass

    def test_2(self):
        reference_coverage = ProbeCoverage(min_depth=131,
                                           percent_coverage=95.2381,
                                           median_depth=155,
                                           k_count=131)
        alt1 = ProbeCoverage(min_depth=1,
                             percent_coverage=100,
                             median_depth=1,
                             k_count=1)
        alternate_coverages = [alt1]
        v1 = VariantProbeCoverage(var_name="A123T",
                                  reference_coverages=[reference_coverage],
                                  alternate_coverages=alternate_coverages
                                  )

        call = self.vt_no_contaim.type(v1)
        self.assertEqual([0, 0], call['genotype'])

    def test_3(self):
        reference_coverage = ProbeCoverage(min_depth=2,
                                           percent_coverage=59.52,
                                           median_depth=2,
                                           k_count=60)
        alt1 = ProbeCoverage(min_depth=1,
                             percent_coverage=83.33,
                             median_depth=1,
                             k_count=83)
        alternate_coverages = [alt1]
        v1 = VariantProbeCoverage(var_name="A123T",
                                  reference_coverages=[reference_coverage],
                                  alternate_coverages=alternate_coverages
                                  )

        call = self.vt2_no_contaim.type(v1)
        self.assertEqual([1, 1], call['genotype'])
        self.assertTrue(call["info"]["conf"] < 100)

    def test_4(self):
        vt = VariantTyper(
            expected_depths=[6],
            contamination_depths=[],
            confidence_threshold=3)
        reference_coverage = ProbeCoverage(min_depth=1,
                                           percent_coverage=100,
                                           median_depth=2,
                                           k_count=2)
        alt1 = ProbeCoverage(min_depth=1,
                             percent_coverage=100,
                             median_depth=1,
                             k_count=1)
        alternate_coverages = [alt1]
        v1 = VariantProbeCoverage(var_name="A123T",
                                  reference_coverages=[reference_coverage],
                                  alternate_coverages=alternate_coverages
                                  )

        call = vt.type(v1)
        self.assertEqual([0, 1], call['genotype'])
        print(call["info"])
        self.assertTrue(call["info"]["conf"] < 100)

import unittest
from atom import Nucleus

class TestNucleus(unittest.TestCase):

    def setUp(self):
        self.u238 = Nucleus('U-238', 238, 92)
        self.pu239 = Nucleus('Pu-239', 239, 94)
        self.cf252 = Nucleus('Cf-252', 252, 98)
        self.te135 = Nucleus('Te-135', 135, 52)

    def test_binding_energy(self):
        """Проверяем расчет энергии связи."""
        self.assertAlmostEqual(self.u238.binding_energy(), 179.88, places=2)
        self.assertAlmostEqual(self.pu239.binding_energy(), 178.72, places=2)
        self.assertAlmostEqual(self.cf252.binding_energy(), 202.80, places=2)

    def test_mass(self):
        """Проверяем расчет массы атома."""
        self.assertAlmostEqual(self.u238.mass(), 238.02891, places=5)
        self.assertAlmostEqual(self.pu239.mass(), 239.05293, places=5)
        self.assertAlmostEqual(self.cf252.mass(), 251.07958, places=5)

    def test_radius(self):
        """Проверяем расчет радиуса атома."""
        self.assertAlmostEqual(self.u238.radius(), 7.44, places=2)
        self.assertAlmostEqual(self.pu239.radius(), 7.45, places=2)
        self.assertAlmostEqual(self.cf252.radius(), 7.47, places=2)

    def test_is_stable(self):
        """Проверяем устойчивость изотопа к бета-распаду."""
        self.assertTrue(self.u238.is_stable())
        self.assertTrue(self.pu239.is_stable())
        self.assertFalse(self.cf252.is_stable()) 

    def test_can_split_even_even(self):
        """Проверяем возможность деления на два четно-четных осколка."""
        self.assertTrue(self.u238.can_split_even_even())
        self.assertTrue(self.pu239.can_split_even_even())
        self.assertTrue(self.te135.can_split_even_even())

if __name__ == '__main__':
    unittest.main()

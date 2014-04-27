#include <vector>
#include <cmath>

#include "Config.hpp"

#include "EigenPimpl.hpp"

#include <boost/filesystem.hpp>

#include "GePolCavity.hpp"
#include "PhysicalConstants.hpp"
#include "Symmetry.hpp"

#include "gtestPimpl.hpp"

namespace fs = boost::filesystem;

class GePolCavityC1Test : public ::testing::Test
{
	protected:
		GePolCavity cavity;
		virtual void SetUp()
		{
  	                Eigen::Vector3d	C1( 0.0000000000, 0.0000000000, 0.0000000000); 	
  		        Eigen::Vector3d	O1( 2.1316110791, 0.0000000000, 0.0000000000);
  		        Eigen::Vector3d	O2(-2.1316110791, 0.0000000000, 0.0000000000);
			std::vector<Sphere> spheres;
			double radiusC = (1.70 * 1.20) / convertBohrToAngstrom;
			double radiusO = (1.52 * 1.20) / convertBohrToAngstrom;
			Sphere sph1(C1, radiusC);
			Sphere sph2(O1, radiusO);
			Sphere sph3(O2, radiusO);
			spheres.push_back(sph1);
			spheres.push_back(sph2);
			spheres.push_back(sph3);
			double area = 0.2 / convertBohr2ToAngstrom2;
			double probeRadius = 1.385 / convertBohrToAngstrom;
			double minRadius = 0.2 / convertBohrToAngstrom;
			// C1
			Symmetry pGroup = buildGroup(0, 0, 0, 0);
			cavity = GePolCavity(spheres, area, probeRadius, minRadius, pGroup);
			fs::rename("PEDRA.OUT", "PEDRA.co2.c1");
			fs::rename("cavity.off", "cavity.co2.c1");
		}
};

TEST_F(GePolCavityC1Test, size)
{
	int size = 448;
	int actualSize = cavity.size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityC1Test, irreducible_size)
{
	int size = 448;
	int actualSize = cavity.irreducible_size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityC1Test, area)
{
	double area = 250.68176442433020;
 	double actualArea = cavity.elementArea().sum();
	EXPECT_NEAR(area, actualArea, 1.0e-10);
}

TEST_F(GePolCavityC1Test, volume)
{
	double volume = 352.55869984340751;
	Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
	Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
	double actualVolume = 0;
        for ( int i = 0; i < cavity.size(); ++i )
	{
		actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
	}
	actualVolume /= 3;
	EXPECT_NEAR(volume, actualVolume, 1.0e-10);
}

class GePolCavityCsTest : public ::testing::Test
{
	protected:
		GePolCavity cavity;
		virtual void SetUp()
		{
  	                Eigen::Vector3d	C1( 0.0000000000, 0.0000000000, 0.0000000000); 	
  		        Eigen::Vector3d	O1( 2.1316110791, 0.0000000000, 0.0000000000);
  		        Eigen::Vector3d	O2(-2.1316110791, 0.0000000000, 0.0000000000);
			std::vector<Sphere> spheres;
			double radiusC = (1.70 * 1.20) / convertBohrToAngstrom;
			double radiusO = (1.52 * 1.20) / convertBohrToAngstrom;
			Sphere sph1(C1, radiusC);
			Sphere sph2(O1, radiusO);
			Sphere sph3(O2, radiusO);
			spheres.push_back(sph1);
			spheres.push_back(sph2);
			spheres.push_back(sph3);
			double area = 0.2 / convertBohr2ToAngstrom2;
			double probeRadius = 1.385 / convertBohrToAngstrom;
			double minRadius = 0.2 / convertBohrToAngstrom;
			// Cs as generated by Oyz
			Symmetry pGroup = buildGroup(1, 1, 0, 0);
			cavity = GePolCavity(spheres, area, probeRadius, minRadius, pGroup);
			fs::rename("PEDRA.OUT", "PEDRA.co2.cs");
			fs::rename("cavity.off", "cavity.co2.cs");
		}
};

TEST_F(GePolCavityCsTest, size)
{
	int size = 448;
	int actualSize = cavity.size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityCsTest, irreducible_size)
{
	int size = 224;
	int actualSize = cavity.irreducible_size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityCsTest, area)
{
	double area = 250.68176442433020;
 	double actualArea = cavity.elementArea().sum();
	EXPECT_NEAR(area, actualArea, 1.0e-10);
}

TEST_F(GePolCavityCsTest, volume)
{
	double volume = 352.55869984340751;
	Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
	Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
	double actualVolume = 0;
        for ( int i = 0; i < cavity.size(); ++i )
	{
		actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
	}
	actualVolume /= 3;
	EXPECT_NEAR(volume, actualVolume, 1.0e-10);
}

class GePolCavityC2Test : public ::testing::Test
{
	protected:
		GePolCavity cavity;
		virtual void SetUp()
		{
  	                Eigen::Vector3d	C1( 0.0000000000, 0.0000000000, 0.0000000000); 	
  		        Eigen::Vector3d	O1( 2.1316110791, 0.0000000000, 0.0000000000);
  		        Eigen::Vector3d	O2(-2.1316110791, 0.0000000000, 0.0000000000);
			std::vector<Sphere> spheres;
			double radiusC = (1.70 * 1.20) / convertBohrToAngstrom;
			double radiusO = (1.52 * 1.20) / convertBohrToAngstrom;
			Sphere sph1(C1, radiusC);
			Sphere sph2(O1, radiusO);
			Sphere sph3(O2, radiusO);
			spheres.push_back(sph1);
			spheres.push_back(sph2);
			spheres.push_back(sph3);
			double area = 0.2 / convertBohr2ToAngstrom2;
			double probeRadius = 1.385 / convertBohrToAngstrom;
			double minRadius = 0.2 / convertBohrToAngstrom;
			// C2 as generated by C2z 
			Symmetry pGroup = buildGroup(1, 3, 0, 0);
			cavity = GePolCavity(spheres, area, probeRadius, minRadius, pGroup);
			fs::rename("PEDRA.OUT", "PEDRA.co2.c2");
			fs::rename("cavity.off", "cavity.co2.c2");
		}
};

TEST_F(GePolCavityC2Test, size)
{
	int size = 448;
	int actualSize = cavity.size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityC2Test, irreducible_size)
{
	int size = 224;
	int actualSize = cavity.irreducible_size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityC2Test, area)
{
	double area = 250.68176442433020;
 	double actualArea = cavity.elementArea().sum();
	EXPECT_NEAR(area, actualArea, 1.0e-10);
}

TEST_F(GePolCavityC2Test, volume)
{
	double volume = 352.55869984340751;
	Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
	Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
	double actualVolume = 0;
        for ( int i = 0; i < cavity.size(); ++i )
	{
		actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
	}
	actualVolume /= 3;
	EXPECT_NEAR(volume, actualVolume, 1.0e-10);
}

class GePolCavityCiTest : public ::testing::Test
{
	protected:
		GePolCavity cavity;
		virtual void SetUp()
		{
  	                Eigen::Vector3d	C1( 0.0000000000, 0.0000000000, 0.0000000000); 	
  		        Eigen::Vector3d	O1( 2.1316110791, 0.0000000000, 0.0000000000);
  		        Eigen::Vector3d	O2(-2.1316110791, 0.0000000000, 0.0000000000);
			std::vector<Sphere> spheres;
			double radiusC = (1.70 * 1.20) / convertBohrToAngstrom;
			double radiusO = (1.52 * 1.20) / convertBohrToAngstrom;
			Sphere sph1(C1, radiusC);
			Sphere sph2(O1, radiusO);
			Sphere sph3(O2, radiusO);
			spheres.push_back(sph1);
			spheres.push_back(sph2);
			spheres.push_back(sph3);
			double area = 0.2 / convertBohr2ToAngstrom2;
			double probeRadius = 1.385 / convertBohrToAngstrom;
			double minRadius = 0.2 / convertBohrToAngstrom;
			// Ci as generated by i
			Symmetry pGroup = buildGroup(1, 7, 0, 0);
			cavity = GePolCavity(spheres, area, probeRadius, minRadius, pGroup);
			fs::rename("PEDRA.OUT", "PEDRA.co2.ci");
			fs::rename("cavity.off", "cavity.co2.ci");
		}
};

TEST_F(GePolCavityCiTest, size)
{
	int size = 448;
	int actualSize = cavity.size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityCiTest, irreducible_size)
{
	int size = 224;
	int actualSize = cavity.irreducible_size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityCiTest, area)
{
	double area = 250.68176442433020;
 	double actualArea = cavity.elementArea().sum();
	EXPECT_NEAR(area, actualArea, 1.0e-10);
}

TEST_F(GePolCavityCiTest, volume)
{
	double volume = 352.55869984340751;
	Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
	Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
	double actualVolume = 0;
        for ( int i = 0; i < cavity.size(); ++i )
	{
		actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
	}
	actualVolume /= 3;
	EXPECT_NEAR(volume, actualVolume, 1.0e-10);
}

class GePolCavityC2hTest : public ::testing::Test
{
	protected:
		GePolCavity cavity;
		virtual void SetUp()
		{
  	                Eigen::Vector3d	C1( 0.0000000000, 0.0000000000, 0.0000000000); 	
  		        Eigen::Vector3d	O1( 2.1316110791, 0.0000000000, 0.0000000000);
  		        Eigen::Vector3d	O2(-2.1316110791, 0.0000000000, 0.0000000000);
			std::vector<Sphere> spheres;
			double radiusC = (1.70 * 1.20) / convertBohrToAngstrom;
			double radiusO = (1.52 * 1.20) / convertBohrToAngstrom;
			Sphere sph1(C1, radiusC);
			Sphere sph2(O1, radiusO);
			Sphere sph3(O2, radiusO);
			spheres.push_back(sph1);
			spheres.push_back(sph2);
			spheres.push_back(sph3);
			double area = 0.2 / convertBohr2ToAngstrom2;
			double probeRadius = 1.385 / convertBohrToAngstrom;
			double minRadius = 0.2 / convertBohrToAngstrom;
			// C2h as generated by Oxy and i
			Symmetry pGroup = buildGroup(2, 4, 7, 0);
			cavity = GePolCavity(spheres, area, probeRadius, minRadius, pGroup);
			fs::rename("PEDRA.OUT", "PEDRA.co2.c2h");
			fs::rename("cavity.off", "cavity.co2.c2h");
		}
};

TEST_F(GePolCavityC2hTest, size)
{
	int size = 448;
	int actualSize = cavity.size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityC2hTest, irreducible_size)
{
	int size = 112;
	int actualSize = cavity.irreducible_size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityC2hTest, area)
{
	double area = 250.68176442433020;
 	double actualArea = cavity.elementArea().sum();
	EXPECT_NEAR(area, actualArea, 1.0e-10);
}

TEST_F(GePolCavityC2hTest, volume)
{
	double volume = 352.55869984340751;
	Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
	Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
	double actualVolume = 0;
        for ( int i = 0; i < cavity.size(); ++i )
	{
		actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
	}
	actualVolume /= 3;
	EXPECT_NEAR(volume, actualVolume, 1.0e-10);
}

class GePolCavityD2Test : public ::testing::Test
{
	protected:
		GePolCavity cavity;
		virtual void SetUp()
		{
  	                Eigen::Vector3d	C1( 0.0000000000, 0.0000000000, 0.0000000000); 	
  		        Eigen::Vector3d	O1( 2.1316110791, 0.0000000000, 0.0000000000);
  		        Eigen::Vector3d	O2(-2.1316110791, 0.0000000000, 0.0000000000);
			std::vector<Sphere> spheres;
			double radiusC = (1.70 * 1.20) / convertBohrToAngstrom;
			double radiusO = (1.52 * 1.20) / convertBohrToAngstrom;
			Sphere sph1(C1, radiusC);
			Sphere sph2(O1, radiusO);
			Sphere sph3(O2, radiusO);
			spheres.push_back(sph1);
			spheres.push_back(sph2);
			spheres.push_back(sph3);
			double area = 0.2 / convertBohr2ToAngstrom2;
			double probeRadius = 1.385 / convertBohrToAngstrom;
			double minRadius = 0.2 / convertBohrToAngstrom;
			// D2 as generated by C2z and C2x
			Symmetry pGroup = buildGroup(2, 3, 6, 0);
			cavity = GePolCavity(spheres, area, probeRadius, minRadius, pGroup);
			fs::rename("PEDRA.OUT", "PEDRA.co2.d2");
			fs::rename("cavity.off", "cavity.co2.d2");
		}
};

TEST_F(GePolCavityD2Test, size)
{
	int size = 448;
	int actualSize = cavity.size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityD2Test, irreducible_size)
{
	int size = 112;
	int actualSize = cavity.irreducible_size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityD2Test, area)
{
	double area = 250.68176442433020;
 	double actualArea = cavity.elementArea().sum();
	EXPECT_NEAR(area, actualArea, 1.0e-10);
}

TEST_F(GePolCavityD2Test, volume)
{
	double volume = 352.55869984340751;
	Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
	Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
	double actualVolume = 0;
        for ( int i = 0; i < cavity.size(); ++i )
	{
		actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
	}
	actualVolume /= 3;
	EXPECT_NEAR(volume, actualVolume, 1.0e-10);
}

class GePolCavityC2vTest : public ::testing::Test
{
	protected:
		GePolCavity cavity;
		virtual void SetUp()
		{
  	                Eigen::Vector3d	C1( 0.0000000000, 0.0000000000, 0.0000000000); 	
  		        Eigen::Vector3d	O1( 2.1316110791, 0.0000000000, 0.0000000000);
  		        Eigen::Vector3d	O2(-2.1316110791, 0.0000000000, 0.0000000000);
			std::vector<Sphere> spheres;
			double radiusC = (1.70 * 1.20) / convertBohrToAngstrom;
			double radiusO = (1.52 * 1.20) / convertBohrToAngstrom;
			Sphere sph1(C1, radiusC);
			Sphere sph2(O1, radiusO);
			Sphere sph3(O2, radiusO);
			spheres.push_back(sph1);
			spheres.push_back(sph2);
			spheres.push_back(sph3);
			double area = 0.2 / convertBohr2ToAngstrom2;
			double probeRadius = 1.385 / convertBohrToAngstrom;
			double minRadius = 0.2 / convertBohrToAngstrom;
			// C2v as generated by Oyz and Oxz
			Symmetry pGroup = buildGroup(2, 1, 2, 0);
			cavity = GePolCavity(spheres, area, probeRadius, minRadius, pGroup);
			fs::rename("PEDRA.OUT", "PEDRA.co2.c2v");
			fs::rename("cavity.off", "cavity.co2.c2v");
		}
};

TEST_F(GePolCavityC2vTest, size)
{
	int size = 448;
	int actualSize = cavity.size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityC2vTest, irreducible_size)
{
	int size = 112;
	int actualSize = cavity.irreducible_size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityC2vTest, area)
{
	double area = 250.68176442433020;
 	double actualArea = cavity.elementArea().sum();
	EXPECT_NEAR(area, actualArea, 1.0e-10);
}

TEST_F(GePolCavityC2vTest, volume)
{
	double volume = 352.55869984340751;
	Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
	Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
	double actualVolume = 0;
        for ( int i = 0; i < cavity.size(); ++i )
	{
		actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
	}
	actualVolume /= 3;
	EXPECT_NEAR(volume, actualVolume, 1.0e-10);
}

class GePolCavityD2hTest : public ::testing::Test
{
	protected:
		GePolCavity cavity;
		virtual void SetUp()
		{
  	                Eigen::Vector3d	C1( 0.0000000000, 0.0000000000, 0.0000000000); 	
  		        Eigen::Vector3d	O1( 2.1316110791, 0.0000000000, 0.0000000000);
  		        Eigen::Vector3d	O2(-2.1316110791, 0.0000000000, 0.0000000000);
			std::vector<Sphere> spheres;
			double radiusC = (1.70 * 1.20) / convertBohrToAngstrom;
			double radiusO = (1.52 * 1.20) / convertBohrToAngstrom;
			Sphere sph1(C1, radiusC);
			Sphere sph2(O1, radiusO);
			Sphere sph3(O2, radiusO);
			spheres.push_back(sph1);
			spheres.push_back(sph2);
			spheres.push_back(sph3);
			double area = 0.2 / convertBohr2ToAngstrom2;
			double probeRadius = 1.385 / convertBohrToAngstrom;
			double minRadius = 0.2 / convertBohrToAngstrom;
			// D2h as generated by Oxy, Oxz and Oyz
			Symmetry pGroup = buildGroup(3, 4, 2, 1);
			cavity = GePolCavity(spheres, area, probeRadius, minRadius, pGroup);
			fs::rename("PEDRA.OUT", "PEDRA.co2.d2h");
			fs::rename("cavity.off", "cavity.co2.d2h");
		}
};

TEST_F(GePolCavityD2hTest, size)
{
	int size = 448;
	int actualSize = cavity.size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityD2hTest, irreducible_size)
{
	int size = 56;
	int actualSize = cavity.irreducible_size();
	EXPECT_EQ(size, actualSize);
}

TEST_F(GePolCavityD2hTest, area)
{
	double area = 250.68176442433020;
 	double actualArea = cavity.elementArea().sum();
	EXPECT_NEAR(area, actualArea, 1.0e-10);
}

TEST_F(GePolCavityD2hTest, volume)
{
	double volume = 352.55869984340751;
	Eigen::Matrix3Xd elementCenter = cavity.elementCenter();
	Eigen::Matrix3Xd elementNormal = cavity.elementNormal();
	double actualVolume = 0;
        for ( int i = 0; i < cavity.size(); ++i )
	{
		actualVolume += cavity.elementArea(i) * elementCenter.col(i).dot(elementNormal.col(i));
	}
	actualVolume /= 3;
	EXPECT_NEAR(volume, actualVolume, 1.0e-10);
}

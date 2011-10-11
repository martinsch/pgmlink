#ifndef ENERGY_H
#define ENERGY_H

#include <stdexcept>
#include <traxels.h>
#include <cmath>

namespace Tracking {

class BotAppearance {
 public:
    double operator()( const Traxel& ) const { return 1000; }
};

class BotDisappearance {
 public:
    double operator()( const Traxel& ) const { return 1000; }
};

class BotMove {
 public:
    double operator()( const Traxel&, const Traxel& ) const { return 200; }
};

class BotDivision {
 public:
    double operator()( const Traxel&, const Traxel&, const Traxel& ) const { return 700; }
};



class GeometryDivision2 {
    public:
    GeometryDivision2(double mean_div_dist, double min_angle) : mean_div_dist_(mean_div_dist), min_angle_(min_angle) {};
    double operator()(const Traxel& ancestor,
	const Traxel& child1,
	const Traxel& child2) const;
    private:
    double mean_div_dist_, min_angle_;
};

class KasterDivision2 {
 public:
 KasterDivision2(double weight, double div_cost) : w_(weight), div_cost_(div_cost) {};
    double operator()(const Traxel& ancestor,
	const Traxel& child1,
	const Traxel& child2) const;
 private:
   double w_;
   double div_cost_;
 };

class NegLnCellness {
 public:
 NegLnCellness(double weight) : w_(weight) {}
  double operator()( const Traxel& ) const;
 private:
  double w_;
};

class NegLnOneMinusCellness {
 public:
 NegLnOneMinusCellness(double weight) : w_(weight) {}
  double operator()( const Traxel& ) const;
 private:
  double w_;
};



////
//// Legacy Energy Functors
////

    ///
    /// Interfaces to general energy functors
    ///
    class NullaryEnergy {
	public:
        virtual double operator()(const Traxels& prev,
                                  const Traxels& curr) const = 0;
        virtual ~NullaryEnergy() = 0;
    };

    class UnaryEnergy {
	public:
        virtual double operator()(const Traxel& t,
                                  const Traxels& prev,
                                  const Traxels& curr) const = 0;
        virtual ~UnaryEnergy() = 0;
    };

    class BinaryEnergy {
	public:
        virtual double operator()(const Traxel& t1,
                                  const Traxel& t2,
                                  const Traxels& prev,
                                  const Traxels& curr) const = 0;
        virtual ~BinaryEnergy() = 0;
    };

    class TertiaryEnergy {
	public:
        virtual double operator()(const Traxel& t1,
                                  const Traxel& t2,
                                  const Traxel& t3,
                                  const Traxels& prev,
                                  const Traxels& curr) const = 0;
        virtual ~TertiaryEnergy() = 0;
    };


    ///
    /// Primitive fixed cost energy functor
    ///
    class ConstantEnergy : public NullaryEnergy, public UnaryEnergy, public BinaryEnergy, public TertiaryEnergy {
	public:
        ConstantEnergy( double energy = 0. );

        virtual double operator()(const Traxel& t1,
                                  const Traxel& t2,
                                  const Traxel& t3,
                                  const Traxels& prev,
                                  const Traxels& curr) const;

        virtual double operator()(const Traxel& t1,
                                  const Traxel& t2,
                                  const Traxels& prev,
                                  const Traxels& curr) const;

        virtual double operator()(const Traxel& t,
                                  const Traxels& prev,
                                  const Traxels& curr) const;


        virtual double operator()(const Traxels& prev,
                                  const Traxels& curr) const;

        double theEnergy;
    };


    
    class SquaredDistance : public BinaryEnergy {
	public:
	SquaredDistance(std::string localisation_feature = "com") : loc_feat_(localisation_feature) {};

        virtual double operator()(const Traxel& from,
                          const Traxel& to,
                          const Traxels& prev,
                          const Traxels& curr) const;
	private:
	std::string loc_feat_;
    };


    class KasterDivision : public TertiaryEnergy {
	public:
	KasterDivision(double division_cost) : div_cost_(division_cost) {};

        virtual double operator()(const Traxel& ancestor,
                          const Traxel& child1,
                          const Traxel& child2,
                          const Traxels& prev,
                          const Traxels& curr) const;
	
	private:
	double div_cost_;
    };

    class GeometryDivision : public TertiaryEnergy {
	public:
	            virtual double operator()(const Traxel& ancestor,
                          const Traxel& child1,
                          const Traxel& child2,
                          const Traxels& prev,
                          const Traxels& curr) const;
    };

    ////
    //// Cellness based mlinder-type energies
    ////
    class CellnessDivision : public TertiaryEnergy {
	public:
        CellnessDivision( double diffCellness = 1461, double absCellness = 190 );

        virtual double operator()(const Traxel& ancestor,
                                  const Traxel& child1,
                                  const Traxel& child2,
                                  const Traxels& prev,
                                  const Traxels& curr) const;
    private:
        // weight for the difference of daughter cells' cellness
        double param_diff_c;
        // weight for the absolute cellness of the parent cell
        double param_abs_c;
    };



    class CellnessMove : public BinaryEnergy {
	public:
        CellnessMove( double diffCellness = 1140 );

        virtual double operator()(const Traxel& from,
                                  const Traxel& to,
                                  const Traxels& prev,
                                  const Traxels& curr) const;
    private:
        // weight for the difference of the cells' cellness
        double param_diff_c;
    };



    class CellnessDisappearance : public UnaryEnergy {
	public:
        CellnessDisappearance( double disappWeight = 1000 );

        virtual double operator()(const Traxel& from,
                                  const Traxels& prev,
                                  const Traxels& curr) const;
    private:
        // weight for the absolute cellness
        double param_abs_c;
    };



    class CellnessAppearance : public UnaryEnergy {
	public:
        CellnessAppearance( double appWeight = 1000 );

        virtual double operator()(const Traxel& to,
                                  const Traxels& prev,
                                  const Traxels& curr) const;
    private:
        // weight for the absolute cellness
        double param_abs_c;
    };

} /* Namespace Tracking */

#endif /* ENERGY_H */

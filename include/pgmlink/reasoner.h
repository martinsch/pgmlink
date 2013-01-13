/**
   @file
   @ingroup matching
   @brief interface for solvers of the Hypotheses Graph matching problem
*/

#ifndef REASONER_H
#define REASONER_H

namespace pgmlink {

class HypothesesGraph;

class Reasoner {
    public:
    /**
     * Setup interal model for given graph.
     */
    virtual void formulate( const HypothesesGraph& ) = 0;

    /**
     * Solve internal model. 
     */
    virtual void infer() = 0;

    /**
     * Write out conclusion from solutions to interal model as property maps of a graph.
     *
     * In general, the graph is the same as the graph used for formulation, but not necessarily
     * (think of storing just the annotations in an empty property graph).
     */
    virtual void conclude( HypothesesGraph& ) = 0;
 };

} /* namespace pgmlink */

#endif /* REASONER_H */

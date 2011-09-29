#include "bot/InputOutput.hxx"
#include "bot/HypothesisSpace.hxx"
#include "bot/ObjectFeatureExtractor.hxx"
#include "bot/AverageObject.hxx"
#include "bot/SingletsGenerator.hxx"
#include "bot/MultipletsGenerator.hxx"
#include "bot/CPLEXSolverSystem.hxx"
#include "bot/TrackingPredictor.hxx"
#include "vigra/hdf5impex.hxx"
#include "bot/SolutionCoder.hxx"
#include "bot/TrainingData.hxx"
#include "bot/TrackingTrainer.hxx"
#include "bot/HDF5ReaderWriter.hxx"

using namespace bot;

int main()
{
    std::string filename("dcelliq-sequence-training.h5");
    // load the image sequence
    std::vector<Matrix2D > images, segmentations;
    TrainingData training;
    HDF5ReaderWriter::load(filename, images, segmentations);
    std::cout << "****Loading the images/segmentations****" << std::endl;
    HDF5ReaderWriter::load(filename, training);
    std::cout << "****Loading the training data****" << std::endl;

    // get the context
    Context context(images);
    std::cout << "****Computing the Context****" << std::endl << context << std::endl << std::endl;

    // load the configuration
    HypothesisSpace space("event-configuration-cell.ini");
    EventConfiguration conf = space.configuration();

    // create singlets/muliplets and extract object features
    std::cout << "****Extracting singlets and multiplets****" << std::endl;
    SingletsSequence singlets_vec;
    SingletsSequence avg_singlet_vec;
    MultipletsSequence multiplets_vec;
    SingletsGenerator singletsGenerator;
    MultipletsGenerator multipletsGenerator(conf.k(), conf.d_max());
    ObjectFeatureExtractor extractor(conf.get_feature_names(), context);
    for (int32 indT = 0; indT < static_cast<int32>(images.size()); indT ++) {
        // generate singlets and multiplets
        Singlets singlets = singletsGenerator(images[indT], segmentations[indT]);
        Multiplets multiplets = multipletsGenerator(images[indT], segmentations[indT], singlets);

        // extract features for them
        extractor(singlets);
        extractor(multiplets);
        // save
        singlets_vec.push_back(singlets);
        avg_singlet_vec.push_back(AverageObject::average(singlets));
        multiplets_vec.push_back(multiplets);
        
        std::cout << "#T = " << indT 
            << ": #singlets = " << singlets.size()
            << ": #multiplets = " << multiplets.size() << std::endl;

    }

    // generate hypotheses and extract joint features
    space(singlets_vec, avg_singlet_vec, multiplets_vec);
    const std::vector<FramePair >& framepairs = space.framepairs();

    // parse the training data
    std::cout << "****Parsing the training data****" << std::endl;
    SolutionCoder coder;
    int32 nTr = training.times().size();
    for (int32 ind = 0; ind < nTr; ind ++) {
        int32 time = training.times()[ind];
        std::cout << "****time = " << time << "****" << std::endl;
        const LabelAssociations& association = training.associations()[ind];

        const std::vector<Event >& events = framepairs[time].events();
        const Singlets& singlets1 = singlets_vec[time];
        const Singlets& singlets2 = singlets_vec[time+1];
        const Multiplets& multiplets1 = multiplets_vec[time];
        const Multiplets& multiplets2 = multiplets_vec[time+1];

        Solution solution;
        coder.decode(
            association, 
            events, 
            singlets1, singlets2,
            multiplets1, multiplets2,
            solution);
        training.solutions().push_back(solution);
    }

    // start the training
    TrackingTrainer trainer;
    const std::vector<Matrix2D > null_vector;
    std::vector<Matrix2D > weights = conf.weights(0.5);
    std::string msg = trainer(training, framepairs, weights, true);
    std::cout << "Training returns: " << msg << std::endl;
    conf.weights() = weights;

    // print the final weights
    std::cout << "Learned weights: " << std::endl;
    conf.print();

    // printe intermediate results: weights, epsilons, losses
    std::cout << "Weights: " << std::endl << trainer.weights() << std::endl << std::endl;
    std::cout << "Epsilons: " << std::endl << trainer.epsilons() << std::endl << std::endl;
    std::cout << "Losses: " << std::endl << trainer.losses() << std::endl << std::endl;

    return 0;
}

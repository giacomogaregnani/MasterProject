include_directories("${PROJECT_SOURCE_DIR}/NumericalSolver")
include_directories("${PROJECT_SOURCE_DIR}/ParticleFilter")

add_library(MCMC
        MCMC.hpp
        MCMC.cpp
        Proposals.hpp
        Proposals.cpp
        Posterior.hpp
        SDEPosterior.hpp
        SDEPosterior.cpp
        PFPosterior.hpp
        PFPosterior.cpp
        PFPosteriorHom.hpp
        PFPosteriorHom.cpp
        NoiselessPosterior.hpp
        NoiselessPosterior.cpp
        ModelingError.hpp
        ModelingError.cpp
        )

target_link_libraries(MCMC
        SDESolver
        ParFil
        )
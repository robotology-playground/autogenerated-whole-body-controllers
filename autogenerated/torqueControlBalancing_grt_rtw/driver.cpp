#ifdef WIN32
#error "Windows is not yet supported"
#endif

#include "torqueControlBalancing.h"
#include <BlockFactory/SimulinkCoder/GeneratedCodeWrapper.h>

#include <csignal>
#include <cstdlib>
#include <iostream>

// Unix dependent
#include <unistd.h>

volatile static bool stop = false;

void termination_handler(int /*sig_number*/)
{
    std::cerr << "CTRL-C pressed. Shutting down gracefully." << std::endl;
    stop = true;
}

int main()
{
    struct sigaction new_action, old_action;
    new_action.sa_handler = termination_handler;
    sigemptyset(&new_action.sa_mask);
    sigaddset(&new_action.sa_mask, SIGTERM);
    new_action.sa_flags = 0;
    sigaction(SIGINT, nullptr, &old_action);
    if (old_action.sa_handler != SIG_IGN) {
        sigaction(SIGINT, &new_action, nullptr);
    }

    blockfactory::coder::GeneratedCodeWrapper<torqueControlBalancingModelClass> model;

    if (!model.initialize()) {
        std::cerr << model.getErrors();
        return EXIT_FAILURE;
    }

    while (true) {
        bool ok = model.step();

        if (!ok || stop) {
            std::cerr << model.getErrors();
            break;
        }
    }

    if (!model.terminate()) {
        std::cerr << model.getErrors();
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

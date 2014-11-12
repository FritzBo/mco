//
//  ep_ui_verifier.cpp
//  mco
//
//  Created by Fritz Bökler on 11.08.14.
//
//

#include "ui_verifier.h"

#include <string>

using std::string;

using TCLAP::CmdLine;
using TCLAP::UnlabeledValueArg;
using TCLAP::ArgException;

#include <mco/generic/benson_dual/upper_image_container.h>

using mco::UpperImageReader;
using mco::UpperImageDDVerifier;
using mco::SimpleUpperImageContainer;

void UpperImageVerifier::perform(int argc, char** argv) {
    try {
        CmdLine cmd("Tool to verify an upper image of an MOLP.", ' ', "0.1");

        UnlabeledValueArg<string> file_name_argument("filename", "Name of the instance file", true, "","filename");

        cmd.add(file_name_argument);

        cmd.parse(argc, argv);

        string file_name = file_name_argument.getValue();

        SimpleUpperImageContainer container;
        UpperImageReader reader;

        reader.read_upper_image(file_name, container);

        UpperImageDDVerifier verifier;

        verifier.verify_double_description(container);

        
    } catch(ArgException& e) {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
}
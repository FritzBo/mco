//
//  ep_ui_verifier.h
//  mco
//
//  Created by Fritz Bökler on 11.08.14.
//
//

#ifndef __mco__ep_ui_verifier__
#define __mco__ep_ui_verifier__

#include "../basic/modules.h"

class UpperImageVerifier : public BasicModule {
    
public:
    virtual void perform(int argc, char** args);
    virtual ~UpperImageVerifier() {}

private:

};

#endif /* defined(__mco__ep_ui_verifier__) */

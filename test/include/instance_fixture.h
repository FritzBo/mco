//
//  instance_fixture.h
//  mco
//
//  Created by Fritz Bökler on 07.08.14.
//
//

#ifndef mco_instance_fixture_h
#define mco_instance_fixture_h

class ParetoInstanceTestFixture
: public ::testing::TestWithParam<tuple<string, unsigned>> {
protected:
    const string filename_ = get<0>(GetParam());
    const unsigned expected_number_of_minimizers_ = get<1>(GetParam());
    
};

#endif

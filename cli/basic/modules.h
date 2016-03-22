//
//  modules.h
//  mco
//
//  Created by Fritz Bökler on 14.07.14.
//
//

#ifndef __mco__modules__
#define __mco__modules__

#include <ostream>
#include <string>
#include <map>
#include <list>
#include <ostream>

#include <tclap/CmdLine.h>

#include <mco/basic/point.h>

enum ModuleType
{
    BASIC,
    ALGORITHM
};

class BasicModule {
public:
    BasicModule()
    :   type(BASIC) {}

    BasicModule(ModuleType mod_type)
    :   type(mod_type) { }

    virtual void perform(int argc, char** args) = 0;
    virtual ~BasicModule() = default;
    const ModuleType type;
};

enum GraphType {
    OGDF,
    FORWARDSTAR
};

class BasicAlgorithmModule : public BasicModule
{
public:
    BasicAlgorithmModule(GraphType type)
    :   BasicModule(ALGORITHM),
        graph_type(type) { }

    const GraphType graph_type;
};

template<class T>
class AlgorithmModule : public BasicAlgorithmModule {
public:
    
    using solution_type = T;
    using csolution_type = const T;
    using solution_type_pointer = T*;
    using csolution_type_pointer = const T*;

    AlgorithmModule(GraphType type)
    :   BasicAlgorithmModule(type) { }
    
    virtual const std::list<std::pair<csolution_type, const mco::Point>>& solutions() = 0;
    virtual std::string statistics() = 0;
};

class ModuleFactory {
public:
    void add_module(std::string name, BasicModule& module);
    
    std::list<std::pair<int, BasicModule*>> parse_module_list(int argc, char** argv);
    
    ~ModuleFactory();

private:
    std::map<std::string, BasicModule*> modules_;
};
#endif /* defined(__mco__modules__) */

#pragma once
/*
 * mcap_parser.h
 *
 *  Created on: 11.11.2013
 *      Author: fritz
 */

#ifndef MCAP_PARSER_H_
#define MCAP_PARSER_H_

#include <string>
#include <memory>

#include <mco/ap/assignment_instance.h>

namespace mco {

class MCAPParser {
public:
	MCAPParser() = delete;

	MCAPParser(std::string filename) :
		filename_(filename) {}

	std::shared_ptr<AssignmentInstance> get_instance();

private:
	std::string filename_;
};

} /* namespace mco */
#endif /* MCAP_PARSER_H_ */

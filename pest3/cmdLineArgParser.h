// -*-c++-*-
// $Id$

// standard includes
#include <string>
#include <cstring>
#include <cstdlib>
#include <map>
#include <vector>

#ifndef CMD_LINE_ARG_PARSER_H
#define CMD_LINE_ARG_PARSER_H

/**
 * @brief Simple command line argument parser
 */
class CmdLineArgParser {

 public:

  /**
   * Constructor
   */
  CmdLineArgParser(int argc, char* argv[]) {

    size_t i = 1;
    while (i < (size_t) argc) {

      std::string sarg(argv[i]);

      // new option?
      if (sarg.size() >= 2 && sarg[0] == '-' 
          && !isdigit(sarg[1]) && sarg[1] != '.') {
        // must be an option
        std::string optionName = sarg.substr(0, 2);
        std::string optionValue = "true"; // options that don't take values
        if (sarg.size() > 2) {
          // no space between name and value
          optionValue = sarg.substr(2);
        }
        else {
          // value is expected to be the next arg
          if (i + 1 < (size_t) argc) {
            std::string sarg2(argv[i + 1]);
            bool nextArgIsOption = true;
            nextArgIsOption &= (sarg2[0] == '-');
            nextArgIsOption &= (!isdigit(sarg2[1]));
            nextArgIsOption &= (sarg2[1] != '.');
            if (!nextArgIsOption) {
              optionValue = sarg2;
              i++;
            }
          }
        }
        std::pair<std::string, std::string> 
          p(optionName, optionValue);
        this->options.insert(p);
      }
      i++;
    }

    this->begin();
  }
  
  /**
   * Destructor
   */
  virtual ~CmdLineArgParser() {}

  /**
   * Set iterator to beginning
   */
  void begin() {
    this->optionsIt = this->options.begin();
  }

  /** 
   * Increment the iterator
   */
  void next() {
    this->optionsIt++;
  }
  
  /**
   * Is iterator valid?
   * @return true if valid false otherwise
   */
  bool ok() {
    return (this->optionsIt != this->options.end());
  }

  // Accessors
  
  /**
   * Get curent option name
   * @return name
   */
  std::string getOptionName() const {
    return this->optionsIt->first;
  }
  
  /**
   * Get curent option value
   * @return name
   */
  std::string getOptionValue() const {
    return this->optionsIt->second;
  }

  /**
   * Retrieve option value 
   * @param name option name
   * @return value, or empty string if invalid name
   * @note: will return the string "true" for options 
   * were set but take no value.
   */
  std::string get(const std::string& name) const {
    std::map<std::string, std::string>::const_iterator it;
    it = this->options.find(name);
    if (it != this->options.end()) {
      return it->second;
    }
    else {
      return std::string("");
    }
  }

  /**
   * Convert string value to list of integers
   * @param name option name
   * @return array
   */
  std::vector<int> getVectInt(const std::string& name) const {
    std::string val = this->get(name);
    if (val.size() > 0) {
      std::vector<int> res;
      int ielem;
      std::string elem;
      bool foundSomething = false;
      size_t i = 0;
      while (i < val.size()) {
        if (foundSomething) {
          elem.push_back(val[i]);
          if (val[i] == ' ') {
            ielem = atoi(elem.c_str());
            res.push_back(ielem);
            elem.resize(0);
            foundSomething = false;
          }
        }
        else {
          if (val[i] != ' ') {
            // detected something
            foundSomething = true;
            elem.push_back(val[i]);
          }
        }
        i++;
      }
      if (elem.size() > 0) {
        ielem = atoi(elem.c_str());
        res.push_back(ielem);
      }
      return res;
    }
    else {
      // invalid name
      return std::vector<int>();
    }
  }
  
 private:
  
  /** options container */
  std::map<std::string, std::string> options;

  /** options iterator */
  std::map<std::string, std::string>::const_iterator optionsIt;

};

#endif // CMD_LINE_ARG_PARSER_H

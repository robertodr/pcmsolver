#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <ctime>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>

#include "BuildInfo.hpp"
#include "LoggerImpl.hpp"

namespace logging
{
    /*! \brief Returns date and time
     */                                      
    inline std::string getTime() {
        std::string time_str;
        time_t raw_time;
                                          
        std::time(&raw_time);
        time_str = std::ctime(&raw_time);
                                          
        // Without the newline character
        return time_str;
    }

    template<typename logPolicy>
    class logger
    {
    private:
	int globalPrintLevel_;
        std::stringstream logStream_;
        logPolicy * policy_;
        std::mutex writeMutex_;
	std::string logFile_;
	bool initialized_;
    
        /*! @name Core printing functionality
         *
         *  A variadic template is used, we specify the version
         *  with a variable number of parameters/types and the version
         *  with no parameters/types.
         *  The variadic template is called recursively.
         *  The type for the first parameter is resolved and streamed
         *  to logStream_. When all the parameters have been streamed
         *  the version with no arguments is called.
         */
        /// @{
        /*! */
        void printImpl() {
            policy_->write(logStream_.str());
            logStream_.str("");
        }

        template<typename First, typename...Rest>
        void printImpl(First parm1, Rest...parm) {
	    logStream_.precision(std::numeric_limits<double>::digits10);
            logStream_ << parm1;
            printImpl(parm...);
        }
        /// @}
       /*! Constructor
        *  \param[in] print the print level
	*  
	*  The build parameters are logged first
         */
        logger(int print = coarse) 
		: globalPrintLevel_(print), policy_(new logPolicy), initialized_(false) 
	{
            if(!policy_) {
                throw std::runtime_error("LOGGER: Unable to create the logger instance");
            }
	    // Write the logfile header
	    logStream_ << "\t\tPCMSolver execution log\n" 
		       << buildInfo() << "\n\t\tLog started : " << getTime() << std::endl;
        }

        /// Destructor
        ~logger() {
            if(policy_) {
                policy_->close_ostream();
                delete policy_;
            }
        }


    public:
        /// User interface for the logger class
        template<int printLvl, typename...Args>
        void print(Args...args) {
	    if (!initialized_) {
	       throw std::runtime_error("Logger not initialized!");
	    }
	    if (globalPrintLevel_ >= printLvl) {
               writeMutex_.lock();   
               printImpl(args...);
               writeMutex_.unlock();
	    }
        }
        static logger& Instance() {
            static logger<FileLogPolicy> loggerInstance;
            return loggerInstance;
        }
	/// Call this to initialize the logger
	void initialize(const std::string & fname = "", int printLvl = 1) {
	    if (!initialized_) {
		if (fname.empty()) {
                    std::stringstream namestream;                                 	
            	    srand(time(NULL));
            	    namestream << "pcmsolver" << "_" << rand() << "_" << getpid();
		    logFile_ = namestream.str();
		} else {
		    logFile_ = fname;
		}
                policy_->open_ostream(logFile_);
		globalPrintLevel_ = printLvl;
		initialized_ = true;
	    }
	}
    };
} // close namespace logging

#endif // LOGGER_HPP

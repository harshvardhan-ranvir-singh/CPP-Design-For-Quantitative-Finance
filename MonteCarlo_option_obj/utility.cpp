//
// Created by Harshvardhan Singh on 12/06/2021.
//


//	utility.cpp


#include "utility.h"
#include <iostream>
#include <stdexcept>



//	PauseAndReturn()
int UtilityFunctions::PauseAndReturn()
{
    char quit = '\0';
    while (quit != 'q')
    {
        std::cout << "Enter q to quit: ";
        std::cin >> quit;
    }
    return 0;
}



//	Putter
void UtilityFunctions::OutputLine(const std::string & text, const double & num)
{
    std::cout << text << " " << num << std::endl;
}

void UtilityFunctions::OutputLine(const std::string & text)
{
    std::cout << text << std::endl;
}





//	end


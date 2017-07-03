#ifndef PARSE_COMMAND_LINE_H
#define PARSE_COMMAND_LINE_H

#include <iostream>
#include <tuple>
using namespace std;

#include "External/MiniCommander.h"
#include "Globals.h"

void summary();

// parse the command line
void parseCommandLine(int argc, const char* argv[])
{
	MiniCommander cmd(argc, argv, false);

	OptionGroup required(Policy::required, "Required Parameters");
	required.addOption("-f", "Path to data file", "--file");
	cmd.addOptionGroup(required);

	OptionGroup optionals(Policy::optional, "Optional Parameters");
	optionals.addOption("-t", "Threshold", "--threshold");
	optionals.addOption("-a", "Algorithm to apply: A* Search (0) or Exhaustive Search (1)", "--algorithm");
	optionals.addOption("-d", "Maximum dimension to be computed", "--dimension");
	optionals.addOption("-h", "Show info and usage", "--help");
	cmd.addOptionGroup(optionals);


	if (!cmd.checkFlags() || cmd.optionExists("--help") || cmd.optionExists("-h"))
	{
		cmd.printHelpMessage("USAGE:");
		exit(EXIT_FAILURE);
	}

	Globals::inputFileName = cmd.getParameter("-f") + cmd.getParameter("--file");
	if (Globals::inputFileName.empty())
	{
		cerr << "Error: please specify required data file." << endl;
		cmd.printHelpMessage("USAGE:");
		exit(EXIT_FAILURE);
	}

	if (cmd.optionExists("-t") || cmd.optionExists("--threshold"))
	{
		std::string temp_threshold = cmd.getParameter("-t") + cmd.getParameter("--threshold");
		if (temp_threshold.empty())
		{
			cerr << "Error: please specify the threshold." << endl;
			cmd.printHelpMessage("USAGE:");
			exit(EXIT_FAILURE);
		}
		Globals::reduction_threshold = stod(temp_threshold);
		Globals::use_optimal_alg = true;
	}

	if (cmd.optionExists("-a") || cmd.optionExists("--algorithm"))
	{
		std::string temp_alg = cmd.getParameter("-a") + cmd.getParameter("--algorithm");
		if (temp_alg.empty())
		{
			cerr << "Error: please specify which algorithm to apply." << endl;
			cmd.printHelpMessage("USAGE:");
			exit(EXIT_FAILURE);
		}
		Globals::which_alg = stoi(temp_alg);
		Globals::use_optimal_alg = true;

		if (Globals::which_alg >= 2)
		{
			cout << "The algorithm index should be 0 (A* Search) or 1 (Exhaustive Search)." << endl;
			exit(EXIT_FAILURE);
		}
	}

	if (cmd.optionExists("-d") || cmd.optionExists("--dimension"))
	{
		std::string temp_dim = cmd.getParameter("-d") + cmd.getParameter("--dimension");
		if (temp_dim.empty())
		{
			cerr << "Error: please specify the maximum dimension to be computed." << endl;
			cmd.printHelpMessage("USAGE:");
			exit(EXIT_FAILURE);
		}
		Globals::max_dim = stoi(temp_dim);
	}

	summary();
}

void summary()
{
	cout << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	cout << "Input data file:  " << Globals::inputFileName << endl;

	cout << "Use use optimal cycle algorithm:  ";
	if (Globals::use_optimal_alg == false)
		cout << "No" << endl << endl;
	else
	{
		cout << "Yes" << endl;
		
		cout << "Algorithm to apply:  ";
		switch (Globals::which_alg)
		{
		case Globals::Algorithm::HEURISTIC_BASED_ALG:
			cout << "Heuristic-based Algorithm" << endl;
			break;
		case Globals::Algorithm::CLASSICAL_ALG:
			cout << "Classical Annotation Algorithm (Exhaustive Search)" << endl;
			break;
		default:
			cout << "Unexpected parameter for algorithm selection!" << endl;
			exit(EXIT_FAILURE);
			break;
		}

		cout << "Threshold:  " << Globals::reduction_threshold << endl;
	}
	cout << "+++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;
}

#endif // !PARSE_COMMAND_LINE_H


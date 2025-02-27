#include "plotter.hpp"

#include <iostream>
#include <fstream>
#include <iomanip>

namespace NPlotter {

TPlotter::TPlotter(const std::string& filename)
    : Width_(800)
    , Height_(600)
    , ImageName_(filename + ".png")
    , DataFile_(filename + ".csv")
    , CommandFile_(filename + ".plt")
    , DataSize_(std::nullopt)
{
}

void TPlotter::SetXValues(const std::vector<int>& values) {
    assert(!DataSize_.has_value() || DataSize_ == values.size());
    if (!DataSize_.has_value()) {
        DataSize_ = values.size();
    }
    XValues_ = values;
}

void TPlotter::Plot() const {
    SaveDataToFile();
    SaveGnuplotCommandsToFile();
    try {
        std::string command = "gnuplot " + CommandFile_; // Запускаем gnuplot с .plt файлом
        int ret = std::system(command.c_str()); // Используем system() вместо popen()
        if (ret != 0) {
            throw std::runtime_error("Failed to execute gnuplot.");
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    } catch (...) {
        std::cerr << "Unknown error while plotting" << std::endl;
    }
}

void TPlotter::SaveDataToFile() const {
    std::ofstream outFile(DataFile_);
    if (!outFile.is_open()) {
        throw std::ios_base::failure("Failed to open the output file: " + DataFile_);
    }

    for (int i = 0; i < DataSize_; ++i) {
        outFile << std::setw(12) << XValues_[i] << "\t";
        for (const auto& chart : Charts_) {
            outFile << std::setw(12) << chart.GetData()[i] << "\t";
        }
        outFile << "\n";
    }
    outFile.flush();
    outFile.close();
}

std::string TPlotter::GenerateGnuplotCommands() const {
    std::ostringstream commands;
    commands << "set term png size " << Width_ << "," << Height_ << "\n";
    commands << "set output '" << ImageName_ << "'\n";
    commands << "set xrange [" << XValues_.front() << ":" << XValues_.back() << "]\n";
    commands << "plot ";

    int index = 2;
    for (const auto& data : Charts_) {
        commands
            << "'" << DataFile_ << "' using 1:"
            << index++ << " with linespoints title '"
            << data.GetTitle() << "', ";
    }
    std::string result = commands.str();
    // remove extra comma
    result.pop_back();
    result.pop_back();
    return result;
}

void TPlotter::SaveGnuplotCommandsToFile() const {
    std::ofstream outFile(CommandFile_);
    if (!outFile.is_open()) {
        throw std::ios_base::failure("Failed to open the Gnuplot command file: " + CommandFile_);
    }

    std::string commands = GenerateGnuplotCommands();
    outFile << commands;
    outFile.close();
}

} // namespace NPlotter

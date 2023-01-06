#ifndef DATAFRAME
#define DATAFRAME
#include "DataFrame.h"
#include <string>
#include <vector>
#include <fstream>

// Implementing Sample
Sample::Sample(double mean) {
    set_mean(mean);
    set_std(0.);
    set_num_samples(1);
}

Sample::Sample(double mean, double std, int num_samples) {
    set_mean(mean);
    set_std(std);
    set_num_samples(num_samples);
}

double Sample::get_mean() {
    return this->mean;
}

void Sample::set_mean(double mean) {
    this->mean = mean;
}

double Sample::get_std() {
    return this->std;
}

void Sample::set_std(double std) {
    this->std = std;
}

double Sample::get_num_samples() {
    return this->num_samples;
}

void Sample::set_num_samples(int num_samples) {
    this->num_samples = num_samples;
}


Sample Sample::combine(Sample* other) {
    int combined_samples = this->num_samples + other->get_num_samples();
    double combined_mean = (this->get_num_samples()*this->get_mean() + other->get_num_samples()*other->get_mean())/combined_samples;
    double combined_std = (((this->get_num_samples()*this->get_std()).pow(2) + (this->get_mean() - combined_mean).pow(2))
                         + ((other->get_num_samples()*other->get_std()).pow(2) + (other->get_mean() - combined_mean).pow(2))
                          )/combined_samples;

    return Sample(combined_mean, combined_std, combined_samples);
}


// Implementing DataSlide

DataSlide::DataSlide() {
    this->data_int = std::map<std::string, int>();
    this->data_double = std::map<std::string, double>();
    this->data = std::map<std::string, std::vector<double>>();
}

void DataSlide::add_int(std::string s, int i) {
    this->data_int.emplace(s, i);
}

void DataSlide::add_double(std::string s, double f) {
    this->data_double.emplace(s, f);
}

void DataSlide::add_data(std::string s) {
    this->data.emplace(s, std::vector<double>(0));
}

void DataSlide::push_data(std::string s, Sample sample) {
    this->data[s].push_back(sample);
}

// Implementing DataFrame

DataFrame::DataFrame(std::vector<DataSlide> slides) {
    for (int i = 0; i < slides.size(); i++) {
        add_slide(slides[i]);
    }
}

DataFrame::DataFrame() {
    this->slides = std::vector<DataSlide>();
}

void DataFrame::add_slide(DataSlide ds) {
    this->slides.push_back(ds);
}

void DataFrame::save(std::string filename) {
    std::ofstream output_file(filename);

    for (int i = 0; i < this->slides.size(); i++) {
        DataSlide slide = this->slides[i];
        for (int j = 0; j < slide.data_int.size(); j++) {
            output_file << 
        }
    }


    output_file.close();
}

}


#endif
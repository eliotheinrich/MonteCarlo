#include "DataFrame.h"
#include <fstream>
#include <iostream>
#include <cmath>

// Implementing Sample
Sample::Sample() {}

Sample::Sample(double mean) : mean(mean), std(0.), num_samples(1) {}

Sample::Sample(double mean, double std, int num_samples) : mean(mean), std(std), num_samples(num_samples) {}


double Sample::get_mean() const {
    return this->mean;
}

void Sample::set_mean(double mean) {
    this->mean = mean;
}

double Sample::get_std() const {
    return this->std;
}

void Sample::set_std(double std) {
    this->std = std;
}

double Sample::get_num_samples() const {
    return this->num_samples;
}

void Sample::set_num_samples(int num_samples) {
    this->num_samples = num_samples;
}


Sample Sample::combine(Sample* other) const {
    int combined_samples = this->num_samples + other->get_num_samples();
    double combined_mean = (this->get_num_samples()*this->get_mean() + other->get_num_samples()*other->get_mean())/combined_samples;
    double combined_std = ((std::pow(this->get_num_samples()*this->get_std(), 2) + std::pow(this->get_mean() - combined_mean, 2))
                         + (std::pow(other->get_num_samples()*other->get_std(), 2) + std::pow(other->get_mean() - combined_mean, 2))
                          )/combined_samples;

    return Sample(combined_mean, combined_std, combined_samples);
}

std::string Sample::to_string() {
    std::string s = "[";
    s += std::to_string(this->mean) + " " + std::to_string(this->std) + " " + std::to_string(this->num_samples) + "]";
    return s;
}

// Implementing DataSlide

DataSlide::DataSlide() {
    this->data_int = std::map<std::string, int>();
    this->data_double = std::map<std::string, double>();
    this->data = std::map<std::string, std::vector<Sample>>();
}

bool DataSlide::contains_int(std::string s) {
    return this->data_int.count(s);
}

bool DataSlide::contains_double(std::string s) {
    return this->data_double.count(s);
}

bool DataSlide::contains_data(std::string s) {
    return this->data.count(s);
}

int DataSlide::get_int(std::string s) {
    if (this->contains_int(s)) {
        return this->data_int[s];
    } else { // TODO better error handling
        return -1;
    }
}

double DataSlide::get_double(std::string s) {
    if (this->contains_double(s)) {
        return this->data_double[s];
    } else {
        return -1.;
    }
}

std::vector<Sample>* DataSlide::get_data(std::string s) {
    if (this->contains_data(s)) {
        return &this->data[s];
    } else {
        return nullptr;
    }
}

void DataSlide::add_int(std::string s, int i) {
    this->data_int.emplace(s, i);
}

void DataSlide::add_double(std::string s, double f) {
    this->data_double.emplace(s, f);
}

void DataSlide::add_data(std::string s) {
    this->data.emplace(s, std::vector<Sample>(0));
}

void DataSlide::push_data(std::string s, Sample sample) {
    this->data[s].push_back(sample);
}

std::string DataSlide::to_string() {
    std::string s = "";
    
    std::map<std::string, int>::iterator int_it = this->data_int.begin();
    s += "\t\t\t" + int_it->first + ": {\"Int\": " + std::to_string(int_it->second) + "}";
    while (this->data_int.end() != ++int_it) {
        s += ",\n\t\t\t" + int_it->first + ": {\"Int\": " + std::to_string(int_it->second) + "}";
    }

    if (!this->data_double.empty()) { s += ",\n"; }
    std::map<std::string, double>::iterator double_it = this->data_double.begin();
    s += "\t\t\t" + double_it->first + ": {\"Float\": " + std::to_string(double_it->second) + "}";
    while (this->data_double.end() != ++double_it) {
        s += ",\n\t\t\t" + double_it->first + ": {\"Float\": " + std::to_string(double_it->second) + "}";
    }

    if (!this->data.empty()) { s += ",\n"; }
    std::map<std::string, std::vector<Sample>>::iterator data_it = this->data.begin();
    s += "\t\t\t" + data_it->first + ": {\"Data\": [";
    std::vector<Sample>::iterator sample_it = data_it->second.begin();
    s += sample_it->to_string();
    while (data_it->second.end() != ++sample_it) {
        s += ", " + sample_it->to_string();
    }
    s += "]}";

    while (this->data.end() != ++data_it) {
        s += ",\n\t\t\t" + data_it->first + ": {\"Data\": [";
        sample_it = data_it->second.begin();
        s += sample_it->to_string();
        while (data_it->second.end() != ++sample_it) {
            s += "," + sample_it->to_string();
        }
    }

    s += "]}";


    return s;
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

void DataFrame::add_param(std::string s, int i) {
    this->params.emplace(s, i);
}

void DataFrame::save(std::string filename) {
    std::string s = "";

    s += "{\n\t\"params\": {\n";

    std::map<std::string, int>::iterator it = this->params.begin();
    s += "\t\t" + it->first + ": {\"Int\": " + std::to_string(it->second) + "}"; 
    while (this->params.end() != ++it) {
        s += ",\n\t\t" + it->first + ": {\"Int\": " + std::to_string(it->second) + "}";
    }

    s += "\n\t},\n\t\"slides\": [\n";

    std::vector<DataSlide>::iterator slides_it = this->slides.begin();
    s += "\t\t{\n" + (*slides_it).to_string() + "\n\t\t}";
    while (this->slides.end() != ++slides_it) {
        s += ",\n\t\t{\n" + (*slides_it).to_string() + "\n\t\t}";
    }

    s += "\n\t]\n}";
    std::cout << s << std::endl;

    // Save to file
    std::ofstream output_file(filename);
    output_file << s;
    output_file.close();
}
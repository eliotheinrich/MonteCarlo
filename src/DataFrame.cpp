#include "DataFrame.h"
#include <fstream>
#include <iostream>
#include <cmath>

std::map<std::string, Sample> to_sample(std::map<std::string, double> *s) {
    std::map<std::string, Sample> new_s;
    for (std::map<std::string, double>::iterator it = s->begin(); it != s->end(); ++it) {
        new_s.emplace(it->first, Sample(it->second));
    }
    return new_s;
};

// Implementing Sample
Sample::Sample() : mean(0.), std(0.), num_samples(0) {}

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
    double combined_std = std::pow(((this->get_num_samples()*std::pow(this->get_std(), 2) + std::pow(this->get_mean() - combined_mean, 2))
                         + (other->get_num_samples()*std::pow(other->get_std(), 2) + std::pow(other->get_mean() - combined_mean, 2))
                          )/combined_samples, 0.5);

    return Sample(combined_mean, combined_std, combined_samples);
}

std::string Sample::to_string() const {
    std::string s = "[";
    s += std::to_string(this->mean) + ", " + std::to_string(this->std) + ", " + std::to_string(this->num_samples) + "]";
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
    this->data.emplace(s, std::vector<Sample>());
}

void DataSlide::push_data(std::string s, Sample sample) {
    this->data[s].push_back(sample);
}

std::string join(const std::vector<std::string> &v, const std::string &delim) {
    std::string s = "";
    for (const auto& i : v) {
        if (&i != &v[0]) {
            s += delim;
        }
        s += i;
    }
    return s;
}

std::string DataSlide::to_string() const {
    std::string s = "\t\t\t";
    

    std::string delim = ",\n\t\t\t";
    std::vector<std::string> buffer;

    for (const auto &[key, val] : this->data_int) {
        buffer.push_back("\"" + key + "\": {\"Int\": " + std::to_string(val) + "}");
    }

    s += join(buffer, delim);

    buffer.clear();
    if (!this->data_double.empty()) {
        s += ",\n\t\t\t";
    }

    for (const auto &[key, val] : this->data_double) {
        buffer.push_back("\"" + key + "\": {\"Float\": " + std::to_string(val) + "}");
    }

    s += join(buffer, delim);

    buffer.clear();
    if (!this->data.empty()) {
        s += ",\n\t\t\t";
    }

    std::vector<std::string> sample_buffer;
    std::string ss;

    for (const auto &[key, samples] : this->data) {
        ss = "\"" + key + "\": {\"Data\": [";
        sample_buffer.clear();
        for (const Sample sample : samples) {
            sample_buffer.push_back(sample.to_string());
        }
        ss += join(sample_buffer, ", ") + "]}";

        buffer.push_back(ss);
    }

    s += join(buffer, delim);

    return s;
}

// Implementing DataFrame

DataFrame::DataFrame(std::vector<DataSlide> slides) {
    this->params = std::map<std::string, int>();
    for (int i = 0; i < slides.size(); i++) {
        add_slide(slides[i]);
    }
}

DataFrame::DataFrame() {
    this->params = std::map<std::string, int>();
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

    std::vector<std::string> buffer;
    for (auto const &[key, val] : this->params) {
        buffer.push_back("\t\t\"" + key + "\": {\"Int\": " + std::to_string(val) + "}");
    }

    s += join(buffer, ",\n");

    s += "\n\t},\n\t\"slides\": [\n";

    buffer.clear();
    int num_slides = this->slides.size();
    for (int i = 0; i < num_slides; i++) {
        buffer.push_back("\t\t{\"data\": {\n" + this->slides[i].to_string() + "}\n\t\t}");
    }

    s += join(buffer, ",\n");

    s += "\n\t]\n}\n";

    // Save to file
    std::ofstream output_file(filename);
    output_file << s;
    output_file.close();
}
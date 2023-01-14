#ifndef DATAFRAME_H
#define DATAFRAME_H

#include <string>
#include <map>
#include <vector>

class Sample {
    private:
        double mean;
        double std;
        double num_samples;

	public:
		Sample();
        Sample(double mean);
		Sample(double mean, double std, int num_samples);

        double get_mean() const;
		void set_mean(double mean);

        double get_std() const;
		void set_std(double std);

        double get_num_samples() const;
		void set_num_samples(int num_samples);

        Sample combine(Sample* other) const;
		std::string to_string();
};

class DataSlide {
	public:
		std::map<std::string, int> data_int;
		std::map<std::string, double> data_double;
		std::map<std::string, std::vector<Sample>> data;

		DataSlide();
		bool contains_int(std::string s);
		bool contains_double(std::string s);
		bool contains_data(std::string s);

		int get_int(std::string s);
		double get_double(std::string s);
		std::vector<Sample>* get_data(std::string s); 
		
		void add_int(std::string s, int i);
		void add_double(std::string s, double f);
		void add_data(std::string s);
		void push_data(std::string s, Sample sample);
		std::string to_string();
};

class DataFrame {
	private:
		std::map<std::string, int> params;
		std::vector<DataSlide> slides;


	public:
		DataFrame();
		DataFrame(std::vector<DataSlide> slides);
		void add_slide(DataSlide ds);
		void add_param(std::string s, int i);
		void save(std::string filename);
};


#endif
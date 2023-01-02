
class DataSlide {
	private:
		std::map<std::string, int> data_int;
		std::map<std::string, double> data_double;
		std::map<std::string, std::vector<double>> data;

	
	public:
		DataSlide();
		void add_int(std::string s, int i);
		void add_double(std::string s, double f);
		void add_data(std::string s);
		void push_data(std::string s, double f);
		
}

class DataFrame {
	private:
		std::vector<DataSlide> slides;


	public:
		DataFrame();
		void add_slide(DataSlide ds);
		void save(std::string filename);
}
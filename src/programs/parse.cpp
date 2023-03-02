namespace gedmap_parse{
	using namespace std;
	using namespace sdsl;
	int main(int argc, char** argv){

		sys_timer t;
		t.start();

		using namespace gedmap_io;
		print_prog_headline("GEDMAP PARSE");
		string fname_fa, fname_vcf, fname_geds;

		gedmap_parse::handle_input(argc, argv, fname_fa, fname_vcf, fname_geds);


		//OUTPUT INFO
		print_row("Parse FA: "	, argv[2]);
		print_row("and  VCF: "	, argv[3]);
		print_row("to  GEDS: "	, fname_geds);
		dotline();

		GEDS_builder gedsb(fname_geds);
		gedsb.build(fname_fa,fname_vcf);

		//WRITE INFO
		dotline();
		print_row("TIME in total: ", t.stop_and_get() , " s" );
		dotline();

		return 0;
	}
}

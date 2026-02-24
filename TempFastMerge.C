void TempFastMerge(const char* inputPattern, const char* outputFile) {
	// Create the TChain for the tree named "tree"
	TChain *chain = new TChain("tree");

	// Add files
	std::cout << "Adding files from: " << inputPattern << std::endl;
	int nFiles = chain->Add(inputPattern);

	if (nFiles == 0) {
		std::cout << "ERROR: No files found! Check your input path." << std::endl;
		return;
	}
	std::cout << "Total files added: " << nFiles << std::endl;

	// Merge using "fast" option (no unzipping)
	std::cout << "Merging to: " << outputFile << " ..." << std::endl;
	chain->Merge(outputFile, "fast");

	std::cout << "Merge complete." << std::endl;
}

#ifndef __TRAITS_HPP__
#define __TRAITS_HPP__

namespace MyTraits{

	//Eigen Types
	//double
	using MatRow        = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
	using MatCol        = Eigen::MatrixXd;
	using VecCol        = Eigen::VectorXd;
	using VecRow        = Eigen::RowVectorXd;
	//unsigned int
	using MatUnsRow     = Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
	using MatUnsCol     = Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
	using VecUnsCol     = Eigen::Matrix<unsigned int, Eigen::Dynamic, 1>;
	using VecUnsRow     = Eigen::Matrix<unsigned int, 1, Eigen::Dynamic>;
	//int
	using MatIntRow     = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
	using MatIntCol     = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
	using VecIntCol     = Eigen::Matrix<int, Eigen::Dynamic, 1>;
	using VecIntRow     = Eigen::Matrix<int, 1, Eigen::Dynamic>;

	//Iterators
	using vector_d_iterator 	= std::vector<double>::iterator;
	using vector_d_citerator 	= std::vector<double>::const_iterator;
	using vector_i_iterator 	= std::vector<int>::iterator;
	using vector_i_citerator 	= std::vector<int>::const_iterator;
	using vector_ui_iterator 	= std::vector<unsigned int>::iterator;
	using vector_ui_citerator 	= std::vector<unsigned int>::const_iterator;


}

#endif

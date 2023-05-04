// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Phil Bradley

// Unit Headers
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/TrialCounter.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh> //Pretty output.

// Numeric Headers
#include <numeric/numeric.functions.hh>
#include <numeric/random/random.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <basic/prof.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/mc.OptionKeys.gen.hh>

#include <protocols/moves/MonteCarloExceptionConverge.hh>
#include <utility/vector1.hh>

// for rosetta++ like boinc graphics
#ifdef BOINC_GRAPHICS
#include <protocols/boinc/boinc.hh>
#endif

///......................................................................
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/simple_moves/SymmetricFragmentMover.hh>
#include <protocols/simple_moves/GunnCost.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/WhileMover.hh>
#include <protocols/abinitio/AllResiduesChanged.hh>
#include <protocols/moves/CompositionMover.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/ozstream.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/templates.OptionKeys.gen.hh>
#include <cstdlib>
#include <string>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/import_pose/import_pose.hh>
#include <utility/io/izstream.hh>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <vector>
#include <iomanip>
#include <string>

using namespace std;
using core::Real;
core::pose::Pose native_pose;
core::scoring::ScoreFunctionOP SCORE3;
int serial_number;
int sequence_lenth;
int Spicker_All_num;
int Spicker_All_file_num;
int Spicker_Number;
ofstream SPICKER_All_Data;
ofstream processdata_1;
ofstream processdata_2;
int N;
int Cluster_number;
bool Save_pdb;
std::map<string, string> parametersMap;
std::string Cluster_score;
///......................................................................

using basic::Error;
using basic::Warning;
static basic::Tracer TR("protocols.moves.MonteCarlo");

namespace protocols
{
	namespace moves
	{

		using namespace core;
		using namespace ObjexxFCL::format;

		/// @details The copy constructor does not copy the OPs, but rather creates new objects using the copy
		/// constructors or the clone() methods of the objects being pointed at.  This is important, since otherwise,
		/// a copy of a Monte Carlo object could corrupt the state held in the original MonteCarlo object.
		MonteCarlo::MonteCarlo(MonteCarlo const &src) : utility::pointer::ReferenceCount(),
														last_accepted_pose_(src.last_accepted_pose_ ? new core::pose::Pose(*src.last_accepted_pose_) : nullptr),
														lowest_score_pose_(src.lowest_score_pose_ ? new core::pose::Pose(*src.lowest_score_pose_) : nullptr),
														temperature_(src.temperature_),
														score_function_(src.score_function_ ? src.score_function_->clone() : core::scoring::ScoreFunctionOP(nullptr)),
														autotemp_(src.autotemp_),
														quench_temp_(src.quench_temp_),
														last_accept_(src.last_accept_),
														mc_accepted_(src.mc_accepted_),
														counter_(src.counter_),
														update_boinc_(src.update_boinc_),
														total_score_of_last_considered_pose_(src.total_score_of_last_considered_pose_),
														last_accepted_score_(src.last_accepted_score_),
														last_score_(src.last_score_),
														lowest_score_(src.lowest_score_),
														heat_after_cycles_(src.heat_after_cycles_),
														convergence_checks_(src.convergence_checks_),
														last_check_(src.last_check_),
														check_frequency_(src.check_frequency_)
		{
		}

		// constructor for monte_carlo object
		MonteCarlo::MonteCarlo(
			Pose const &init_pose,		   // PoseCOP init_pose,
			ScoreFunction const &scorefxn, // ScoreFunctionCOP scorefxn,
			Real const temperature) : temperature_(temperature),
									  autotemp_(false),
									  quench_temp_(0.0),
									  last_accept_(0),
									  mc_accepted_(MCA_accepted_score_beat_last), // init_pose beats the absence of a pose
									  counter_(TrialCounterOP(new TrialCounter)),
									  update_boinc_(true),
									  total_score_of_last_considered_pose_(0.0),
									  last_accepted_score_(0.0),
									  last_score_(0.0),
									  lowest_score_(0.0),
									  heat_after_cycles_(150)
		{
			last_accepted_pose_ = PoseOP(new Pose());
			lowest_score_pose_ = PoseOP(new Pose());
			// score_function_ = new ScoreFunction(scorefxn);
			score_function_ = scorefxn.clone();
			reset(init_pose);

			last_check_ = 0;
			check_frequency_ = basic::options::option[basic::options::OptionKeys::mc::convergence_check_frequency]();
		}

		MonteCarlo::MonteCarlo(
			ScoreFunction const &scorefxn, // ScoreFunctionCOP scorefxn,
			Real const temperature) : temperature_(temperature),
									  autotemp_(false),
									  quench_temp_(0.0),
									  last_accept_(0),
									  mc_accepted_(MCA_accepted_score_beat_last), // an empty pose beats the absence of a pose
									  counter_(TrialCounterOP(new TrialCounter)),
									  update_boinc_(true),
									  total_score_of_last_considered_pose_(0.0),
									  last_accepted_score_(0.0),
									  last_score_(0.0),
									  lowest_score_(0.0),
									  heat_after_cycles_(150)
		{
			last_accepted_pose_ = PoseOP(new Pose());
			lowest_score_pose_ = PoseOP(new Pose());
			// score_function_ = new ScoreFunction(scorefxn);
			score_function_ = scorefxn.clone();
			last_check_ = 0;
			check_frequency_ = basic::options::option[basic::options::OptionKeys::mc::convergence_check_frequency]();
		}

		MonteCarlo::~MonteCarlo() = default;

		void MonteCarlo::clear_poses()
		{
			last_accepted_pose_ = PoseOP(new Pose());
			lowest_score_pose_ = PoseOP(new Pose());
		}

		void
		MonteCarlo::reset_scorefxn(
			Pose const &init_pose,
			ScoreFunction const &scorefxn)
		{
			// score_function_ = new ScoreFunction(scorefxn);
			score_function_ = scorefxn.clone();
			reset(init_pose);
		}

		void
		MonteCarlo::set_temperature(Real const temp)
		{
			temperature_ = temp;
		}

		/// return the simulation state to the lowest energy structure we've seen
		void
		MonteCarlo::recover_low(Pose &pose)
		{
			(pose) = (*lowest_score_pose_);
			*last_accepted_pose_ = *lowest_score_pose_;
			// last_accepted_pose_ = new Pose( *lowest_score_pose_ );
			last_accepted_score_ = last_accepted_pose_->energies().total_energy();
		}

		void
		MonteCarlo::show_state() const
		{
			TR << "MC: " << temperature_
			   << "  " << (*score_function_)(*last_accepted_pose_)
			   << "  " << (*score_function_)(*lowest_score_pose_)
			   << "  " << last_accepted_score()
			   << "  " << lowest_score()
			   << "  " << last_accept_
			   << "  " << autotemp_
			   << "  " << quench_temp_
			   << "  " << to_string(mc_accepted_)
			   << std::endl;
			show_counters();
		}

		/////////////////////////////////////////////////////////////////////////////

		void
		MonteCarlo::show_scores() const
		{
			TR << "MonteCarlo:: last_accepted_score,lowest_score: " << last_accepted_score() << ' ' << lowest_score() << std::endl;
		}
		/////////////////////////////////////////////////////////////////////////////
		void
		MonteCarlo::reset_counters()
		{
			counter_->reset();
		}

		/// @detail return number of trials since last reset
		core::Size
		MonteCarlo::total_trials() const
		{
			return counter_->total_trials();
		}
		/////////////////////////////////////////////////////////////////////////////
		void
		MonteCarlo::show_counters() const
		{
			counter_->show();
		}

		void
		MonteCarlo::change_weight(core::scoring::ScoreType const &t, Real const &setting)
		{
			score_function_->set_weight(t, setting);
		}

		/// set the scorefxn,  re-scores last-accepted and lowest-score pose
		void
		MonteCarlo::score_function(ScoreFunction const &scorefxn)
		{
			using namespace scoring;

			// *score_function_ = scorefxn;
			// score_function_ = new ScoreFunction(scorefxn);
			score_function_ = scorefxn.clone(); // JAB - this is dangerous if the scorefunction has not been allocated on the heap.
			/// Why not just accept a ScoreFunctionCOP like almost every other class in Rosetta???????

			// TR << "new score_function_ within mc!" << std::endl;

			lowest_score_ = (*score_function_)(*lowest_score_pose_);
			// TR << "lowest_score: " << lowest_score_ << " total: " << lowest_score_pose_->energies().total_energy() << std::endl;
			/// Now handled automatically.  score_function_->accumulate_residue_total_energies( *lowest_score_pose_ );

			if (false)
			{ // DEBUG
				pose::Pose copy_pose;
				copy_pose = *lowest_score_pose_;
				copy_pose.energies().clear();
				TR << "score copy" << std::endl;
				Real const copy_score = (*score_function_)(copy_pose);
				if (std::abs(copy_score - lowest_score_) > 1E-6)
				{
					TR << "Score discrepancy.  lowest_score: " << lowest_score_ << " vs copy score: " << copy_score << std::endl;
					TR << "pose score: ";
					lowest_score_pose_->energies().total_energies().show_if_nonzero_weight(TR, score_function_->weights());
					TR << std::endl;

					TR << "copy score: ";
					copy_pose.energies().total_energies().show_if_nonzero_weight(TR, score_function_->weights());
					TR << std::endl;

					TR << "Difference: ";
					core::scoring::EnergyMap emap = lowest_score_pose_->energies().total_energies();
					emap -= copy_pose.energies().total_energies();
					emap.show_if_nonzero_weight(TR, score_function_->weights());
					TR << std::endl;
				}
			}

			// T("protocols.moves.MonteCarlo.score_function") << "lowest_score";
			// score_function_->show( T("protocols.moves.MonteCarlo.score_function"), *lowest_score_pose_ );
			// T("protocols.moves.MonteCarlo.score_function");

			last_accepted_score_ = (*score_function_)(*last_accepted_pose_);
			// TR << "last_accepted_score: " << last_accepted_score << " total: " << last_accepted_pose_->energies().total_energy() << std::endl;
			/// Now handled automatically.  score_function_->accumulate_residue_total_energies( *last_accepted_pose_ );

			if (false)
			{ // DEBUG
				pose::Pose copy_pose;
				copy_pose = *last_accepted_pose_;
				copy_pose.energies().clear();
				TR << "score copy" << std::endl;
				Real const copy_score = (*score_function_)(copy_pose);
				if (std::abs(copy_score - last_accepted_score_) > 1E-6)
				{
					TR << "Score discrepancy.  last_accepted_score: " << last_accepted_score_ << " vs copy score: " << copy_score << std::endl;
					TR << "pose score: ";
					last_accepted_pose_->energies().total_energies().show_if_nonzero_weight(TR, score_function_->weights());
					TR << std::endl;

					TR << "copy score: ";
					copy_pose.energies().total_energies().show_if_nonzero_weight(TR, score_function_->weights());
					TR << std::endl;

					TR << "Difference: ";
					core::scoring::EnergyMap emap = last_accepted_pose_->energies().total_energies();
					emap -= copy_pose.energies().total_energies();
					emap.show_if_nonzero_weight(TR, score_function_->weights());
					TR << std::endl;
				}
			}

			// T("protocols.moves.MonteCarlo.score_function") << "last_accepted";
			// score_function_->show( T("protocols.moves.MonteCarlo.score_function"), *last_accepted_pose_ );
			// T("protocols.moves.MonteCarlo.score_function") << " finished!";

			/// SML 6/25/08 - if last accepted is now better than lowest, change it!
			/// explicit this-> necessary to resolve function name against identically named local variables
			if (this->last_accepted_score() < this->lowest_score())
			{
				PROF_START(basic::MC_ACCEPT);
				*lowest_score_pose_ = *last_accepted_pose_;
				lowest_score_ = last_accepted_score_;
				PROF_STOP(basic::MC_ACCEPT);
			}
		}

		/////////////////////////////////////////////////////////////////////////////
		//////////////////////////////
		// @details behavior depends on setting of temperature
		//
		// return true for an accept, false otherwise
		//
		//  mc_accepted
		//   3 = accepted:score beat low score and last_accepted score
		//   2 = accepted:score beat last_accepted score
		//   1 = thermally accepted: score worse than last_accepted score
		//   0 = not accepted
		//
		// Optional inputs:
		//
		//  proposal_density_ratio = ratio of backward proposal probability to
		//                            forward proposal probability,
		//                            to maintain detailed balance.
		//
		//  inner_score_delta_over_temperature
		//                          = change in energy in any separate inner loop that
		//                            used Boltzmann criterion with different energy
		//                            function. Don't penalize with that energy
		//                            difference again. See, e.g.,
		//                             Hetenyi et al., JCP 117: 8203-8207.
		//                            ( Note: functionally redundant with
		//                               proposal_density_ratio)

		bool
		MonteCarlo::boltzmann(
			Pose &pose,
			std::string const &move_type,						// = unk
			core::Real const proposal_density_ratio,			// = 1
			core::Real const inner_score_delta_over_temperature // = 0
		)
		{

			// Work around a current bug in the pose observer classes..
#ifdef BOINC_GRAPHICS
			if (update_boinc_)
			{
				boinc::Boinc::update_graphics_current(pose);
			}
#endif

			Real score = ((*score_function_)(pose));
			;
			// now delegate decision making...
			bool const accept(boltzmann(score, pose, move_type, proposal_density_ratio, inner_score_delta_over_temperature));

			return accept; // accept!
		}

		//////////////////////////////////////////////////////////
		// See notes above on boltzmann()
		//////////////////////////////////////////////////////////
		bool
		MonteCarlo::boltzmann(
			core::Real score,
			std::string const &move_type,						 // = "unk"
			core::Real const proposal_density_ratio,			 // = 1
			core::Real const inner_score_delta_over_temperature, // = 0
			bool check_lowest_score								 //=true
		)
		{
			// figure out the actual score
			total_score_of_last_considered_pose_ = score; // save this for the TrialMover so that it may keep statistics.
			counter_->count_trial(move_type);

#ifdef BOINC_GRAPHICS
			if (update_boinc_)
			{
				boinc::Boinc::update_mc_trial_info(counter_->trial(move_type), move_type);
			}
#endif

			last_score_ = score;

			Real const score_delta(score - last_accepted_score_);
			Real const boltz_factor_score = (-score_delta / temperature_) + inner_score_delta_over_temperature;
			Real const probability = std::exp(std::min(40.0, std::max(-40.0, boltz_factor_score))) * proposal_density_ratio;

			if (probability < 1)
			{
				if (numeric::random::rg().uniform() >= probability)
				{
					mc_accepted_ = MCA_rejected; // rejected

					autotemp_reject();
					return false; // rejected
				}
				mc_accepted_ = MCA_accepted_thermally; // accepted thermally
			}
			else
			{
				mc_accepted_ = MCA_accepted_score_beat_last; // energy is lower than last_accepted
			}

			counter_->count_accepted(move_type);
			counter_->count_energy_drop(move_type, score_delta);
			last_accepted_score_ = score;

			autotemp_accept();

			if (check_lowest_score && score < lowest_score())
			{
				lowest_score_ = score;
				mc_accepted_ = MCA_accepted_score_beat_low; // 3;
			}
			return true; // accept!
		}

		bool
		MonteCarlo::boltzmann(
			core::Real score,
			Pose &pose,
			std::string const &move_type,						 // = "unk"
			core::Real const proposal_density_ratio,			 // = 1
			core::Real const inner_score_delta_over_temperature, // = 0
			bool check_lowest_score								 //=true
		)
		{
			bool accept = boltzmann(score, move_type, proposal_density_ratio, inner_score_delta_over_temperature, check_lowest_score);

			//============================================================
			// save spicker file and pose information
			Real pscore = (*SCORE3)(pose);

			if (accept)
			{

				Real dm_score = DM_score(native_pose, pose);

				// outpu information
				Output_SPICKER_All_Data(pose, score);
				cout <<"*****************" << std::endl;
				processdata_1
					<< Spicker_Number << '\t' << setw(5) << setprecision(3)
					<< core::scoring::CA_rmsd(native_pose, pose) << '\t'
					<< dm_score << '\t'
					<< pscore << std::endl;

				if (Save_pdb)
				{
					// save .pdb file
					stringstream pose_name;
					pose_name
						<< "./data/" << Spicker_All_file_num
						<< "/pose/pose_" << Spicker_Number << ".pdb";

					pose.dump_pdb(pose_name.str());
				}

				processdata_2 << "pose " << Spicker_Number << std::endl;
				processdata_2 << pscore << std::endl;

				// output angle data
				for (int i = 1; i <= sequence_lenth; ++i)
					processdata_2 << pose.phi(i) << setw(15) << pose.psi(i) << std::endl;
			}

			//=================================================================================

			// rejected ?
			if (!accept)
			{
				evaluate_convergence_checks(pose, true /*reject*/, false /* not final*/);
				pose = (*last_accepted_pose_);
				return false; // rejected
			}

			// accepted !
			PROF_START(basic::MC_ACCEPT);
			*last_accepted_pose_ = pose;

			// print out the scores for each decoy to cmd out, if you pass a flag
			// nice for testing
			if (basic::options::option[basic::options::OptionKeys::mc::log_scores_in_MC].user() && basic::options::option[basic::options::OptionKeys::mc::log_scores_in_MC]() == true)
			{
				score_function_->show(pose);
			}

#ifdef BOINC_GRAPHICS
			if (update_boinc_)
			{
				boinc::Boinc::update_graphics_last_accepted(pose, last_accepted_score());
			}
#endif

			if (mc_accepted_ == MCA_accepted_score_beat_low)
			{
				*lowest_score_pose_ = pose;
				evaluate_convergence_checks(pose, false /*not reject*/, false /*not final*/);

#ifdef BOINC_GRAPHICS
				if (update_boinc_)
				{
					boinc::Boinc::update_graphics_low_energy(pose, lowest_score());
				}
#endif
			}

			PROF_STOP(basic::MC_ACCEPT);

			return true;
		}

		void
		MonteCarlo::reset(Pose const &pose)
		{
			PROF_START(basic::MC_ACCEPT);
			*last_accepted_pose_ = pose;
			PROF_STOP(basic::MC_ACCEPT);

			Real const score = (*score_function_)(*last_accepted_pose_);

			/// Now handled automatically.  score_function_->accumulate_residue_total_energies( *last_accepted_pose_ );

			PROF_START(basic::MC_ACCEPT);
			*lowest_score_pose_ = *last_accepted_pose_;
			PROF_STOP(basic::MC_ACCEPT);

			last_score_ = score;
			last_accepted_score_ = score;
			lowest_score_ = score;
		}

		/////////////////////////////////////////////////////////////////////////////
		void
		MonteCarlo::set_autotemp(
			bool const setting,
			Real const quench_temp)
		{
			autotemp_ = setting;
			quench_temp_ = quench_temp;
			last_accept_ = 0;
		}

		void
		MonteCarlo::set_last_accepted_pose(Pose const &pose)
		{
			*last_accepted_pose_ = pose;
			last_accepted_score_ = last_accepted_pose_->energies().total_energy();
		}

		void MonteCarlo::set_last_accepted_pose(core::pose::Pose const &pose, core::Real score)
		{
			*last_accepted_pose_ = pose;
			last_accepted_score_ = score;
		}

		void MonteCarlo::set_lowest_score_pose(core::pose::Pose const &pose)
		{
			*lowest_score_pose_ = pose;
			lowest_score_ = pose.energies().total_energy();
		}

		void MonteCarlo::set_lowest_score_pose(core::pose::Pose const &pose, core::Real score)
		{
			*lowest_score_pose_ = pose;
			lowest_score_ = score;
		}

		bool
		MonteCarlo::eval_lowest_score_pose(
			Pose &pose,
			bool score_pose,			 // true
			bool update_stats,			 // false
			std::string const &move_type // unk
		)
		{
			// Get or calculate energy
			Real score;
			if (score_pose)
			{
				score = (*score_function_)(pose);
			}
			else
			{
				score = pose.energies().total_energy();
			}

			last_score_ = score;

			// Evaluate
			total_score_of_last_considered_pose_ = score; // save this for the TrialMover so that it may keep statistics.
			if (score < lowest_score())
			{
				*lowest_score_pose_ = pose;
				lowest_score_ = score;
				if (update_stats)
				{
					counter_->count_accepted(move_type);
					counter_->count_energy_drop(move_type, score - last_accepted_score());
					last_accepted_score_ = score;
					mc_accepted_ = MCA_accepted_score_beat_low;
					*last_accepted_pose_ = pose;
					evaluate_convergence_checks(pose, false /*not reject*/, false /*not final*/);
				}

				return true;
			}
			else
			{
				if (update_stats)
				{
					mc_accepted_ = MCA_rejected; // rejected
				}
				return false;
			}
		}

		core::scoring::ScoreFunction const &
		MonteCarlo::score_function() const
		{
			return *score_function_;
		}

		Real
		MonteCarlo::last_accepted_score() const
		{
			// if (last_accepted_pose_->energies().total_energy() != last_accepted_score_) {
			//  TR << "Last: " << last_accepted_pose_->energies().total_energy() << " != " << last_accepted_score_ << " diff " << last_accepted_pose_->energies().total_energy() - last_accepted_score_ << std::endl;
			// }
			// return last_accepted_pose_->energies().total_energy();
			return last_accepted_score_;
		}

		Real
		MonteCarlo::last_score() const
		{
			return last_score_;
		}

		void
		MonteCarlo::set_last_score(core::Real last_score)
		{
			last_score_ = last_score;
		}

		Real
		MonteCarlo::lowest_score() const
		{
			// if (lowest_score_pose_->energies().total_energy() != lowest_score_) {
			//  TR << "Low: " << lowest_score_pose_->energies().total_energy() << " != " << lowest_score_ << " diff " << lowest_score_pose_->energies().total_energy() - lowest_score_ << std::endl;
			// }
			// return lowest_score_pose_->energies().total_energy();
			return lowest_score_;
		}

		MCA
		MonteCarlo::mc_accepted() const
		{
			return mc_accepted_;
		}

		std::string
		MonteCarlo::mc_accepted_string() const
		{
			return to_string(mc_accepted_);
		}

		/////////////////////////////////////////////////////////////////////////////
		// replicate logic from monte_carlo.cc
		//
		// should probably make these parameters ( 150, 1.0 )
		void
		MonteCarlo::autotemp_reject()
		{
			// int const heat_after_cycles( 150 );
			Real const heat_delta(quench_temp_ * 0.5);
			Real const max_temperature(quench_temp_ * 10.0);

			if (!autotemp_)
				return;
			if (last_accept_ >= (int)heat_after_cycles_)
			{
				// if ( temperature_ > max_temperature * 0.25 )
				TR << "autotemp_reject -- heat: " << last_accept_ << ' ' << temperature_ << std::endl;
				last_accept_ = -1;
				set_temperature(std::min(temperature_ + heat_delta, max_temperature));
			}
			++last_accept_;
		}

		/////////////////////////////////////////////////////////////////////////////
		// replicate logic from monte_carlo.cc
		void
		MonteCarlo::autotemp_accept()
		{
			if (!autotemp_)
				return;
			if (temperature_ != quench_temp_)
			{
				set_temperature(quench_temp_);
				TR << "autotemp_accept: reset temperature_ = " << temperature_ << std::endl;
			}

			last_accept_ = 0;
		}

		void
		MonteCarlo::evaluate_convergence_checks(core::pose::Pose const &pose, bool reject, bool /*final*/)
		{
			if (!reject || numeric::mod(last_check_++, check_frequency_) == 0)
			{
				for (auto &convergence_check : convergence_checks_)
				{
					(*convergence_check)(pose, *this, reject);
				}
			}
		}

		void
		MonteCarlo::push_back(moves::MonteCarloExceptionConvergeOP check)
		{
			convergence_checks_.push_back(check);
		}

		// for Python bindings
		std::ostream &operator<<(std::ostream &os, MonteCarlo const &mc)
		{
			os << "protocols.moves.MonteCarlo:\n";
			os << "\"Temperature\" (kT): " << mc.temperature() << "\n";
			os << "Total Trials: " << mc.total_trials() << "\n";
			os << "Lowest Score: " << mc.lowest_score() << "\n";
			os << "Last Accepted Score: " << mc.last_accepted_score() << "\n";
			os << "last_accept = " << mc.lowest_score() << std::endl;
			return os;
		}

		Real MonteCarlo::DM_score(core::pose::Pose &pose, core::pose::Pose &tempPose)
		{
			Real score_sum = 0;
			Real dm_score = 0;
			Real num = 0;
			Real d0 = 0;
			Real varepsilon = 0.001;
			for (int i = 0; i < sequence_lenth;)
			{
				for (int j = i + 2; j < sequence_lenth;)
				{
					Vector a1 = (pose.residue(i + 1).atom("CA").xyz());
					Vector a2 = (pose.residue(j + 1).atom("CA").xyz());

					Vector b1 = (tempPose.residue(i + 1).atom("CA").xyz());
					Vector b2 = (tempPose.residue(j + 1).atom("CA").xyz());

					d0 = log(varepsilon + fabs(i - j));
					Real di = fabs(a1.distance(a2) - b1.distance(b2));
					if (fabs(i - j) >= 3)
					{
						score_sum += 1 / (1 + pow(di / d0, 2));
						num += 1;
					}
					j = j + 2;
				}
				i = i + 2;
			}
			dm_score = score_sum / num;

			return dm_score;
		}

		void
		MonteCarlo::mkdir(string &file_name)
		{
			string path;
			if (Save_pdb == true)
				path = "mkdir -p " + file_name + "/pose";

			else
				path = "mkdir -p " + file_name;

			cout << path << endl;
			system(path.c_str());

			string tempfile;
			tempfile = "cp ./tools/spicker " + file_name;
			system(tempfile.c_str());

			tempfile = "chmod 755 " + file_name;
			tempfile.append("/spicker");
			system(tempfile.c_str());
		}

		void MonteCarlo::Output_SPICKER_All_Data(core::pose::Pose pose, Real energy)
		{
			++Spicker_All_num;
			++Spicker_Number;

			if (Spicker_All_num > 13000)
			{
				Spicker_All_num = 1;
				string file_name;
				stringstream ss;

				SPICKER_Demand_All(native_pose);
				++Spicker_All_file_num;
				ss << "./data/" << Spicker_All_file_num;
				ss >> file_name;
				mkdir(file_name);

				ss.clear();
				ss << "./data/" << Spicker_All_file_num << "/spicker.data1";
				ss >> file_name;
				SPICKER_All_Data.open(file_name.c_str());
			}

			int Num_resdue(pose.total_residue());
			SPICKER_All_Data << right << setw(8) << Num_resdue
							 << right << setw(10) << setprecision(3) << fixed << energy
							 << right << setw(8) << Spicker_All_num
							 << right << setw(8) << Spicker_All_num
							 << endl;
			for (int r = 1; r <= Num_resdue; ++r)
			{
				numeric::xyzVector<Real> CA_ = pose.residue(r).xyz("CA");
				SPICKER_All_Data << right << setw(10) << setprecision(3) << fixed << CA_.x()
								 << right << setw(10) << setprecision(3) << fixed << CA_.y()
								 << right << setw(10) << setprecision(3) << fixed << CA_.z()
								 << endl;
			}
		}

		void MonteCarlo::SPICKER_Demand_All(core::pose::Pose Pose)
		{
			SPICKER_All_Data.close();

			string file_name;
			stringstream ss;

			ss << "./data/" << Spicker_All_file_num << "/tra.in";
			ss >> file_name;
			ofstream tra(file_name.c_str());
			tra << " 1 1 1" << endl;
			tra << "spicker.data1" << endl;
			tra.close();

			ss.clear();
			ss << "./data/" << Spicker_All_file_num << "/rmsinp";
			ss >> file_name;
			ofstream rmsinp(file_name.c_str());
			rmsinp << 1 << "  " << sequence_lenth << "    ! these two numbers indicate region for RMSD calculations [1," << sequence_lenth << "]"
				   << "\n"
				   << sequence_lenth << "       ! protein length" << std::endl;
			rmsinp.close();

			ss.clear();
			ss << "./data/" << Spicker_All_file_num << "/seq.dat";
			ss >> file_name;
			ofstream seq(file_name.c_str());
			for (int r = 1; r <= sequence_lenth; ++r)
				seq << right << setw(5) << r
					<< right << setw(6) << Pose.residue(r).name3()
					<< endl;
			seq.close();

			ss.clear();
			ss << "./data/" << Spicker_All_file_num << "/CA";
			ss >> file_name;
			ofstream CA(file_name.c_str());
			for (int r = 1; r <= sequence_lenth; ++r)
			{
				numeric::xyzVector<Real> CA_ = Pose.residue(r).xyz("CA");

				CA << "ATOM"
				   << right << setw(7) << r
				   << right << setw(4) << "CA"
				   << right << setw(5) << Pose.residue(r).name3()
				   << right << setw(6) << r
				   << right << setw(12) << CA_.x()
				   << right << setw(8) << CA_.y()
				   << right << setw(8) << CA_.z()
				   << endl;
			}
			CA.close();
		}
	} // namespace moves
} // namespace core

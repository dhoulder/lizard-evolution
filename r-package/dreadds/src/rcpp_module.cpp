// -*- mode: C++; indent-tabs-mode: nil; -*-

/**
 * R API for the DREaD_ds model
 */

#include <Rcpp.h>
//[[Rcpp::plugins(cpp11)]]

#include <vector>
#include <string>

#include <model-config.h>
#include <model.h>
#include <species.h>
#include <deme.h>
#include "environment.h"

using namespace Rcpp;

static DreadDs::Model *model_factory(CharacterVector args) {
  std::vector<const char *> av;
  av.push_back("dreadds");
  for (String s: args)
    av.push_back(s.get_cstring());
  DreadDs::Config conf(av.size(), av.data());
  conf.verbosity = 0; // keep std::cerr quiet
  return new DreadDs::Model(conf);
}

static NumericVector get_extent(DreadDs::Model *m) {
  return wrap((m->env).get_extent());
}


/**
 * Return a list of NumericMatrix() representing the environment
 * values at the current step
 */
static List get_env(DreadDs::Model *m) {
  const DreadDs::Config &conf = m->conf;
  const auto &&s = m->env.values.shape();
  DreadDs::EnvCell ec;
  bool no_data;
  NumericMatrix e[conf.env_dims];
  for (int i=0; i < conf.env_dims; ++i)
    e[i] = NumericMatrix(s[0], s[1]);
  DreadDs::Location loc;

  for (loc.y = 0; loc.y < s[0]; ++loc.y)
    for (loc.x = 0; loc.x < s[1]; ++loc.x) {
      m->env.get(loc, ec, &no_data);
      for (int i=0; i < conf.env_dims; ++i) {
        e[i](loc.y, loc.x) = no_data? NA_REAL : ec[i];
      }
    }
  List ret(conf.env_dims);
  for (int i=0; i < conf.env_dims; ++i) {
    ret[i] = e[i];
  }
  return ret;
}


/**
 * Returns demes for all extant species as a list of data
 * frames. Each list element is named with the species name. Each data
 * frame has columns labelled "row", "column", "amount", followed by
 * the environment values for each environment variable, followed by
 * columns describing the genetic position.
 */
static List get_demes(DreadDs::Model *m) {
  bool no_data;
  const DreadDs::Config &conf = m->conf;
  int ncol = 3 + conf.env_dims*3 + conf.genetic_dims;
  CharacterVector col_names(ncol);
  int cni = 0;
  col_names[cni++] = "row";
  col_names[cni++] = "column";
  col_names[cni++] = "amount";
  for (int i=0; i < conf.env_dims; ++i) {
    std::string si = std::to_string(i);
    col_names[cni++] = "env_"+ si;
    col_names[cni++] = "niche_centre_" + si;
    col_names[cni++] = "niche_breadth_" + si;
  }
  for (int i=0; i < conf.genetic_dims; i++)
    col_names[cni++] = "genetic_position_" + std::to_string(i);

  CharacterVector species_names;
  List df_list(m->tips.size());
  int dfli = 0;
  for (auto &&species: m->tips) {
    int n_demes = species->demes->size();

    IntegerVector row(n_demes), column(n_demes);
    NumericVector amount(n_demes);

    NumericVector env[conf.env_dims],
      nc[conf.env_dims],
      nb[conf.env_dims];
    for (int i=0; i < conf.env_dims; ++i) {
      env[i] = NumericVector(n_demes);
      nc[i] = NumericVector(n_demes);
      nb[i] = NumericVector(n_demes);
    }
    NumericVector gp[conf.genetic_dims];
    for (int i=0; i < conf.genetic_dims; i++)
      gp[i] = NumericVector(n_demes);
    species_names.push_back(species->get_name());

    DreadDs::EnvCell ec;
    auto cell_itr = species->demes->begin();
    for (int di = 0; di < n_demes; ++di, ++cell_itr) {
      const DreadDs::Location &loc = cell_itr->first;
      m->env.get(loc, ec, &no_data);
      for (auto &&deme: cell_itr->second) {
	auto &&g = deme.genetics;
        row[di] = loc.y + 1; // Indexes start at 1 in R
        column[di] = loc.x + 1;
        amount[di] = deme.amount;
	for (int i=0; i < conf.env_dims; ++i) {
          env[i][di] = (no_data? NA_REAL : ec[i]);
          nc[i][di] = g.niche_centre[i];
          nb[i][di] = g.niche_tolerance[i] *2.0f; // *2 to convert to breadth
        }
	for (int i=0; i < conf.genetic_dims; i++) {
	  gp[i][di] = g.genetic_position[i];
        }
      }
    }

    // The recommended way to create a DataFrame with a variable
    // number of columns is to build a List() and then convert it. See
    // https://stackoverflow.com/a/8631853 https://stackoverflow.com/a/8648201
    List df_columns(ncol);
    int dfcol = 0;
    df_columns[dfcol++] = row;
    df_columns[dfcol++] = column;
    df_columns[dfcol++] = amount;
    for (int i=0; i < conf.env_dims; ++i) {
      std::string si = std::to_string(i);
      df_columns[dfcol++] = env[i];
      df_columns[dfcol++] = nc[i];
      df_columns[dfcol++] = nb[i];
    }
    for (int i=0; i < conf.genetic_dims; i++) {
      df_columns[dfcol++] = gp[i];
    }
    df_columns.names() = col_names;
    df_list[dfli++] = DataFrame(df_columns);
  }

  df_list.names() = species_names;
  return df_list;
}


/**
 * Return a character vector that describes the phylogeny descended
 * from each initial species using Newick format.  See
 * https://en.wikipedia.org/wiki/Newick_format
 * https://www.rdocumentation.org/packages/ape/versions/5.2/topics/read.tree
 * https://www.r-phylo.org/wiki/HowTo/InputtingData
 * https://www.r-phylo.org/wiki/HowTo/InputtingTrees
 */
static CharacterVector get_phylogeny(DreadDs::Model *m) {
  CharacterVector ret(m->roots.size());
  auto sp_itr = m->roots.begin();
  for (int i=0; i < m->roots.size(); ++i, ++sp_itr) {
    ret[i] = (*sp_itr)->phylogeny_as_newick();
  }
  return ret;
}


/**
 * Return a DataFrame describing the current characteristics of all
 * species. Each row corresponds to a species
 */
static DataFrame get_species(DreadDs::Model *m) {
  const DreadDs::Config &conf = m->conf;
  auto all_species =  m->get_all_species();
  int n = all_species.size();
  IntegerVector ids(n), parents(n), extinctions(n),
    splits(n), steps(n), cell_counts(n);
  NumericVector populations(n);
  CharacterVector species_names(n);

  int n_stats = conf.env_dims * 6;
  int n_gen = conf.genetic_dims * 2;
  NumericVector latest_stats[n_stats];
  NumericVector genetics[n_gen];
  for (int i=0; i < n_stats; ++i)
    latest_stats[i] = NumericVector(n);
  for (int i=0; i < n_gen; i++)
    genetics[i] = NumericVector(n);

  auto sp_itr = all_species.begin();
  for (int sp_i=0; sp_i < all_species.size(); ++sp_i, ++sp_itr) {
    auto &species = *(*sp_itr);
    auto &stats = species.latest_stats;
    ids[sp_i] = species.id;
    species_names[sp_i] = species.get_name();
    parents[sp_i] = species.parent? species.parent->id : -1;
    extinctions[sp_i] = species.extinction;
    splits[sp_i] = species.split;
    steps[sp_i] = species.step;
    cell_counts[sp_i] = stats.cell_count;
    populations[sp_i] = stats.population;

    for (int i=0, col=0; i < conf.env_dims; ++i) {
      latest_stats[col++][sp_i] = stats.niche_stats[i].position_mean;
      latest_stats[col++][sp_i] = stats.niche_stats[i].position_sd;
      latest_stats[col++][sp_i] = stats.niche_stats[i].breadth_mean;
      latest_stats[col++][sp_i] = stats.niche_stats[i].breadth_sd;
      latest_stats[col++][sp_i] = stats.niche_stats[i].max;
      latest_stats[col++][sp_i] = stats.niche_stats[i].min;
    }

    for (int i=0, col=0; i < conf.genetic_dims; i++) {
      genetics[col++][sp_i] = stats.genetic_stats[i].mean;
      genetics[col++][sp_i] = stats.genetic_stats[i].sd;
    }
  }
  int ncol = 7 + n_stats + n_gen;
  List df_columns(ncol);
  CharacterVector col_names(ncol);
  int cni = 0;
  col_names[cni++] = "id";
  col_names[cni++] = "parent";
  col_names[cni++] = "extinction";
  col_names[cni++] = "split";
  col_names[cni++] = "step";
  col_names[cni++] = "cell_count";
  col_names[cni++] = "population";

  for (int i=0; i < conf.env_dims; ++i) {
    std::string si = std::to_string(i);
    col_names[cni++] = "niche_position_mean_" + si;
    col_names[cni++] = "niche_position_sd_" + si;
    col_names[cni++] = "niche_breadth_mean_" + si;
    col_names[cni++] = "niche_breadth_sd_" + si;
    col_names[cni++] = "niche_max_" + si;
    col_names[cni++] = "niche_min_" + si;
  }
  for (int i=0, col=0; i < conf.genetic_dims; i++) {
    std::string si = std::to_string(i);
    col_names[cni++] = "genetic_mean_" + si;
    col_names[cni++] = "genetic_sd_" + si;
  }
  int dfcol = 0;
  df_columns[dfcol++] = ids;
  df_columns[dfcol++] = parents;
  df_columns[dfcol++] = extinctions;
  df_columns[dfcol++] = splits;
  df_columns[dfcol++] = steps;
  df_columns[dfcol++] = cell_counts;
  df_columns[dfcol++] = populations;
  for (int i=0; i < n_stats; ++i)
    df_columns[dfcol++] = latest_stats[i];
  for (int i=0; i < n_gen; i++)
    df_columns[dfcol++] = genetics[i];
  df_columns.names() = col_names;
  DataFrame df(df_columns);
  df.attr("row.names") = species_names;
  return df;
}

/**
 * Run model for a specified number of steps or until it finishes or
 * reached defined iteration limit.
 */
static int run_steps(DreadDs::Model *m, int n) {
  int r = -1;
  for (int i=0; r!=0 && i<n; ++i)
    r = m->do_step();
  return r;
}

static int run_all(DreadDs::Model *m) {
  return run_steps(m, m->conf.n_iterations - m->step);
}


RCPP_MODULE(dreadds){
  using namespace Rcpp ;

  class_<DreadDs::Model>("dreaddsModel")
    // see http://www2.uaem.mx/r-mirror/web/packages/Rcpp/news.html
    // https://github.com/RcppCore/Rcpp/pull/938
    // https://stackoverflow.com/questions/23599043/expose-class-in-rcpp-factory-instead-of-constructor
    .factory<CharacterVector>(model_factory)
    .method("runSteps", &run_steps)
    .method("runAll", &run_all)
    .method("save", &DreadDs::Model::save)
    .method("getEnv", &get_env)
    .method("getExtent", &get_extent)
    .method("getDemes", &get_demes)
    .method("getPhylogeny", &get_phylogeny)
    .method("getSpecies", &get_species)
    .field_readonly("step", &DreadDs::Model::step)
    ;
}

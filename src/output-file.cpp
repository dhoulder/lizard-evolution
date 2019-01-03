// -*- coding: utf-8 -*-

// For C open(), fdopen() etc.
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include <unistd.h>

#include "model-config.h"
#include "exceptions.h"
#include "output-file.h"

/**
 * Using plain old C file IO because C++ doesn't provide an
 * equivalent to O_CREAT|O_EXCL.
 *
 * Note that this assumes that there are no signal handlers, and so
 * doesn't handle errno == EINTR
 */


FILE *DreadDs::open_output_file(const Config &conf,
				std::string filename_suffix) {
  std::string output_filename = conf.output_dir +
    "/" +
    conf.output_file_prefix +
    filename_suffix;

  int fd = open(output_filename.c_str(),
		// atomically create and open, fail if exists
		O_WRONLY | O_CREAT | O_EXCL,
		0664);
  FILE *of = NULL;
  if (fd > -1)
    of = fdopen(fd, "w");
  if (NULL == of) {
    if (fd > -1) close(fd);
    char msg_buf[200] = "";
#if (_POSIX_C_SOURCE >= 200112L) && !  _GNU_SOURCE
    strerror_r(errno, msg_buf, sizeof(msg_buf));
    char *mb = msg_buf;
#else
    char *mb = strerror_r(errno, msg_buf, sizeof(msg_buf));
#endif
    throw ApplicationException("Could not create " +
			       output_filename + ": " + mb);
  }
  return of;
}

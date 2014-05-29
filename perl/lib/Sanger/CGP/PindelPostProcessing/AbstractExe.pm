####################################################
# Copyright (c) 2013 Genome Research Ltd.
# Author: Cancer Genome Project, cgpit@sanger.ac.uk
# See LICENCE.TXT for details
####################################################
package Sanger::CGP::PindelPostProcessing::AbstractExe;

use FindBin;
use Sanger::CGP::Pindel;
use Exporter 'import';
@EXPORT_OK = qw( get_version );

sub get_version {
  return $VERSION;
}

1;

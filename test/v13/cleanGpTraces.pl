use strict;
use warnings;

while (my $line = <>) {
    next if ($line =~ /[Ee]xited with (status)|(exit) code/);  # Local server executor appends this as final line; exact text varies by GP version.

    $line =~ s/\[[\d]{2}:[\d]{2}:[\d]{2}\]/\[timestamp\]/g;
    $line =~ s/(\/{0,1})Axis[0-9]+\.att\_/$1/g;
    $line =~ s/\/[A-Za-z\/0-9@._-]+\//\[system_path\]\//g;
    
    # Clean out R platform details and base install package versions
    $line =~ s/Platform:(.+)\(64-bit\)/Platform:\[R platform\]\(64-bit\)/g;
    $line =~ s/KernSmooth_[-_.0-9]+/KernSmooth/g;
    $line =~ s/codetools_[-_.0-9]+/codetools/g;
    
    print $line;
}

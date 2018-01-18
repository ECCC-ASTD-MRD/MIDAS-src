package htmling;

use strict;

### CONSTANTS
$htmling::dblspace = "  ";
$htmling::indentspace = $htmling::dblspace x 2;
$htmling::headerspace = $htmling::indentspace;
$htmling::comment_indent = $htmling::indentspace x 2;

### PUBLIC GLOBALS
$htmling::comments_type = "smart";
$htmling::suppress_calls = 1;
$htmling::calls_make_links = 0;
$htmling::html_filenames_original_case = 0;

### GLOBALS
$htmling::htmlfile = "";
$htmling::indent = 0;

# Return the name of the HTML file for the specified PROGRAM or MODULE
sub html_filename {
  my ($name) = @_;
  $name = lc $name unless $htmling::html_filenames_original_case;
  return $name . ".html";
}

# This is the main calling point from f90doc.
# Takes all top-level objects: programs, subroutines, functions, and modules.
# Warns if given something else.
sub do_toplevel {
  my ($top, $outfile) = @_;
  
  my $type = $top->{'type'};
  unless ($type eq 'module' || $type eq 'subroutine' || $type eq 'function' ||
          $type eq 'program') {
    warn "Warning: Unrecognized top-level object $type will not be documented.\n";
    return;
  }

  # A positive-length name.  Necessary because programs may not have names.
  if (defined $outfile) {
    $htmling::htmlfile = $outfile . html_filename (
        ($top->{'name'} eq '' ? $type : $top->{'name'}));
  } else {
    $htmling::htmlfile = html_filename (
        ($top->{'name'} eq '' ? $type : $top->{'name'}));
  }
  print "Generating $htmling::htmlfile...\n";
  open OUT, ">$htmling::htmlfile";

  print OUT "<HTML style=\"font-family:arial\">\n";
  print OUT "<HEAD>\n";
  print OUT "   <TITLE> $type $top->{'name'} (generated by f90doc) </TITLE>\n";
  print OUT "</HEAD>\n";
  print OUT "<H1 style=\"font-weight:normal;color:darkred\"> ", lc ($type), " $top->{'name'} </H1>\n";

  do_comments ($top->{'comments'}, 0);
  print OUT "</PRE><HR>";

  print OUT "<PRE STYLE=\"FONT-SIZE:14PX\">$type $top->{'name'}\n";

  list_uses (@{$top->{'uses'}});
  list_calls (1, keys %{$top->{'calls'}}) if exists $top->{'calls'};
#  list_html ("Types", map (($_->{'type'} eq "type" ? ($_) : ()), @{$top->{'ocontains'}}));
  list_html ("Public Variables", map (($_->{'type'} eq "var" ? ($_) : ()), @{$top->{'ocontains'}}));
  list_html ("Interfaces", map (($_->{'type'} eq "interface" ? ($_) : ()), @{$top->{'ocontains'}}));
  list_html ("Subroutines and functions", map (($_->{'type'} eq "subroutine" || $_->{'type'} eq "function" ? ($_) : ()), @{$top->{'ocontains'}}));

  print OUT "\nend $type $top->{'name'}\n";

  #do_comments ($top->{'comments'}, 1);
  print OUT "</PRE>";

  my @list;
#  @list = map (($_->{'type'} eq "type" ? ($_) : ()), @{$top->{'ocontains'}});
#  print OUT "\n<HR><H2> Description of Types </H2>\n" if @list;
#  do_html (@list);
#  @list = map (($_->{'type'} eq "var" ? ($_) : ()), @{$top->{'ocontains'}});
#  print OUT "\n<HR><H2> Description of Variables </H2>\n" if @list;
#  do_html (@list);
#  @list = map (($_->{'type'} eq "interface" ? ($_) : ()), @{$top->{'ocontains'}});
#  print OUT "\n<HR><H2> Description of Interfaces </H2>\n" if @list;
#  do_html (@list);
  @list = map (($_->{'type'} eq "subroutine" || $_->{'type'} eq "function" ? ($_) : ()), @{$top->{'ocontains'}});
  print OUT "\n<HR><H2 style=\"font-weight:normal\"> Description of Subroutines and Functions </H2>\n" if @list;
  do_html (@list);
   
  print OUT "</HTML>\n";
  close OUT;
}

sub list_uses {
  if (@_) {
    print OUT "\n${htmling::indentspace}${htmling::headerspace}! Uses\n";
    my ($use);
    foreach $use (@_) {
      my ($module, $extra) = @$use;
      $extra = defined $extra ? ", $extra" : "";
      print OUT "${htmling::indentspace}",
                "use <A style=\"text-decoration:none\" HREF=\"", html_filename ($module),
                "\">$module</A>$extra\n";
    }
  }
}

sub list_calls {
  return if $htmling::suppress_calls;
  my ($big, @calls) = (@_);
  if (@calls) {
    @calls = sort @calls;
    @calls = map { "<A style=\"text-decoration:none\" HREF=\"$_.html\">$_</A>" } @calls
      if $htmling::calls_make_links;
    if ($big) {
      print OUT join ("\n",
          "\n${htmling::indentspace}${htmling::headerspace}! Calls",
          (map { "${htmling::indentspace}call $_" } @calls), "");
    } else {
      print OUT "${htmling::indentspace}! Calls: ", join (", ", @calls), "\n";
    }
  }
}

sub list_html {
  my ($title) = shift;

  if (@_) {
    print OUT "\n${htmling::indentspace}${htmling::headerspace}! $title\n";
    my ($struct);
    foreach $struct (@_) {
      my ($name, $type) = (txt2html ($struct->{'name'}), $struct->{'type'});
      my ($href) = "<A style=\"text-decoration:none\" HREF=\"${htmling::htmlfile}#${type}_" .
        lc ($name) . "\">$name</A>";
      if ($type eq "var") {
          # only display public variables
	  if ($struct->{'vis'} eq "public") {
              print OUT $htmling::indentspace;
              print OUT var2str ($struct, $href) . "\n";
          }
      } elsif ($type eq "subroutine" ||
               $type eq "function") {
        print OUT $htmling::indentspace;
        print OUT join (" ", attriblist ($struct), "");
        print OUT typing::type_to_f90 ($struct->{'rtype'}) . " "
          if exists $struct->{'rtype'};
        my $flag;
        for $flag ('recursive', 'elemental', 'pure') {
          print OUT "$flag " if $struct->{$flag};
        }
        print OUT "$type $href";
        print OUT " (" . join (", ", @{$struct->{'parms'}}) . ")";
        print OUT " result ($struct->{'result'})"
          if exists $struct->{'result'} && !exists $struct->{'rtype'};
        print OUT "\n";
      } else {
        print OUT $htmling::indentspace;
        print OUT join (" ", attriblist ($struct), "");
        print OUT "$type $href\n";
      }
    }
  }
}

sub do_html {
   if (@_) {
      my ($struct);

      foreach $struct (@_) {
         my ($name, $type) = (txt2html ($struct->{'name'}), $struct->{'type'});
         if (! $htmling::indent) {
            print OUT "<A style=\"text-decoration:none\" NAME=\"${type}_" . lc ($name) .
               "\"><H3 style=\"font-weight:normal;color:darkred\">$name</H3></A>\n";
            print OUT "<PRE STYLE=\"FONT-SIZE:14PX\">";
         }

         print OUT $htmling::indentspace x $htmling::indent;
         if ($type eq "var") {
             print OUT var2str ($struct) . "\n";
         } elsif ($type eq "mprocedure") {
             die "do_html: bare module procedure $struct->{'name'} (no enclosing module)"
                 unless exists $struct->{'bind'};
             print OUT
                 "module procedure <A style=\"text-decoration:none\" HREF=\"#$struct->{'bind'}->{'type'}_" .
                 lc ($struct->{'name'}) . "\">$name</A>\n";
         } elsif ($type eq "subroutine" || $type eq "function") {
             print OUT join (" ", attriblist ($struct), "");
             print OUT typing::type_to_f90 ($struct->{'rtype'}) . " "
                 if exists $struct->{'rtype'} && !exists $struct->{'result'};
             my $flag;
             for $flag ('recursive', 'elemental', 'pure') {
               print OUT "$flag " if $struct->{$flag};
             }
             print OUT "$type $name";
             print OUT " (" . join (", ", @{$struct->{'parms'}}) . ")";
             print OUT " result ($struct->{'result'})"
               if exists $struct->{'result'};
             print OUT "\n";
         } else {
             print OUT join (" ", attriblist ($struct), "");
             print OUT "$type $name\n";
         }

         $htmling::indent++;

         if ($type eq "var" || $type eq "mprocedure") {
         } elsif ($type eq "type") {
           print OUT $htmling::indentspace x $htmling::indent, "private\n"
             if exists $struct->{'privatetype'};
           print OUT $htmling::indentspace x $htmling::indent, "sequence\n"
             if exists $struct->{'sequencetype'};
           do_html (@{$struct->{'ocontains'}});
         } elsif ($type eq "interface") {
           do_html (@{$struct->{'ocontains'}});
         } elsif ($type eq "subroutine" || $type eq "function") {
           my @interest = @{$struct->{'parms'}};
           push @interest, $struct->{'result'} if exists $struct->{'result'};
           push @interest, $name
             if $type eq "function" && !exists $struct->{'result'} &&
               !exists $struct->{'rtype'};
           my $arg;
           foreach $arg (@interest) {
             my (@things) = values %{$struct->{'contains'}->{lc $arg}};
             #die "Confused by/no declaration for parameter $arg of $type $name"
             #  if scalar @things != 1;
             #do_html ($things[0]);
           }
         } else {
           #die "do: I don't know what a $type is";
         }

         list_calls (0, keys %{$struct->{'calls'}}) if exists $struct->{'calls'};

         $htmling::indent--;

         if ($type ne "var" && $type ne "mprocedure") {
            print OUT $htmling::indentspace x $htmling::indent . "end $type $name\n";
         }

         do_comments ($struct->{'comments'}, ! $htmling::indent);
      }
   }
}

# Pass comments and a flag saying if you want to end the current <PRE> block.
sub do_comments {
   my ($comments, $endpre) = @_;
   if ($comments eq "") {
      print OUT "</PRE>\n" if $endpre;
      return;
   }

   #print OUT "\n" unless $htmling::indent;

   if ($htmling::comments_type eq "preformatted") {
      my ($s) = $htmling::indentspace x $htmling::indent . $htmling::comment_indent;
      $comments =~ s/^/$s/m if $htmling::indent;
      $comments =~ s/^\n*//s;
      $comments =~ s/\n*$//s;
      print OUT $comments, "\n";
      print OUT "</PRE>\n" if $endpre;
   } else {
      print OUT "</PRE>\n";
      print OUT "<DL><DD><DL><DD>\n" if $htmling::indent;
      if ($htmling::comments_type eq "html") {
      } elsif ($htmling::comments_type eq "smart") {
         my @newcomments = ();
         my $verbmode = 0;
         my @listmode = ();
         my $line;
         foreach $line (split ("\n", $comments)) {
            if ($verbmode) {
              if ($line =~ /^>/) {
                warn "`$line' found while already in verbatim mode";
                substr ($line, 0, 1) = " ";
                push @newcomments, $line;
              } elsif ($line =~ /^</) {
                $verbmode = 0;
                substr ($line, 0, 1) = " ";
                push @newcomments, $line . "</PRE>";
              } elsif ($line =~ /^v/) {
                warn "`$line' found while already in verbatim mode";
                substr ($line, 0, 1) = " ";
                push @newcomments, $line;
              } else {
                push @newcomments, $line;
              }
              next;
            }

            # _italic_ and *bold*
            while ($line =~ /(\A|\W)_(\w|\w.*?\w)_(\Z|\W)/) {
              my ($left, $mid, $right) = ("$`$1<I>", $2, "</I>$3$'");
              $mid =~ s/_/ /g;
              $line = $left . $mid . $right;
            }
            while ($line =~ /(\A|\W)\*(\w|\w.*?\w)\*(\Z|\W)/) {
              my ($left, $mid, $right) = ("$`$1<STRONG>", $2, "</STRONG>$3$'");
              $mid =~ s/\*/ /g;
              $line = $left . $mid . $right;
            }

            # Lists
            if ($line =~ /^( *)-/) {
              if (! @listmode || length ($1) > $listmode[$#listmode]) {
                push @listmode, length $1;
                push @newcomments, $1 . "<UL>";
              } else {
                while ($listmode[$#listmode] != length ($1)) {
                  push @newcomments, " " x $listmode[$#listmode] . "</UL>";
                  pop @listmode;
                  die "Unindented to invalid position in `$line'"
                    unless @listmode;
                }
              }
              push @newcomments, $1 . "<LI> " . substr ($line, length ($&));
            } elsif ($line =~ /^>/) {
              #warn "Verbatim mode started in list mode" if @listmode;
              $verbmode = 1;
              substr ($line, 0, 1) = " ";
              push @newcomments, "<PRE STYLE=\"FONT-SIZE:14PX\">" . $line;
            # Ignore $line =~ /^</ because it may be an HTML tag.
            } elsif ($line =~ /^v/) {
              #warn "One-line verbatim in list mode" if @listmode;
              substr ($line, 0, 1) = " ";
              push @newcomments, "<PRE STYLE=\"FONT-SIZE:14PX\">$line</PRE>";
            } elsif ($line =~ /^\s*$/) {
              push @newcomments, "<P>";
            } elsif (@listmode) {
              $line =~ /^( *)(\t?)/;
              warn "Tabs have strange effects on indentation detection"
                if length ($2) > 0;
              while (@listmode && $listmode[$#listmode] > length ($1)) {
                push @newcomments, " " x $listmode[$#listmode] . "</UL>";
                pop @listmode;
              }
              push @newcomments, $line;
            } else {
              push @newcomments, $line;
            }
         }
         my $list;
         foreach $list (@listmode) {
             push @newcomments, " " x $list . "</UL>";
         }
         $comments = join ("\n", @newcomments);
      } else {
         die "Unsupported comments type `$htmling::comments_type'";
      }
      $comments =~ s/<P>\n(<P>\n)+/<P>\n/g;
      $comments =~ s/<P>\n$//;
      $comments =~ s/^<P>\n//;
      $comments =~ s/<P>/<DD>/g if $htmling::indent;
      print OUT $comments . "\n";
      print OUT "</DL></DL>\n" if $htmling::indent;
      print OUT "<PRE STYLE=\"FONT-SIZE:14PX\">" unless $endpre;
   }
}

sub var2str {
    my ($var, $href) = @_;

    my ($typestr) = typing::type_to_f90 ($var->{'vartype'});
    my ($initial) = (!exists $var->{'initial'} ? ""
          : " $var->{'initop'} " . typing::expr_to_f90 ($var->{'initial'}));
    $href = txt2html ($var->{'name'}) unless $href;
    return $typestr . join (", ", "", attriblist ($var)) . " :: $href$initial";
}

sub txt2html {
    my ($txt) = @_;
    $txt =~ s/</&lt;/g;
    $txt =~ s/>/&gt;/g;
    return $txt;
}

sub attriblist {
    my ($struct) = @_;
    my @attribs = ();

    push @attribs, $struct->{'vis'} if exists $struct->{'vis'};
    push @attribs, "optional" if exists $struct->{'optional'};
    push @attribs, @{$struct->{'tempattribs'}}
        if exists $struct->{'tempattribs'};

    return @attribs;
}

1;

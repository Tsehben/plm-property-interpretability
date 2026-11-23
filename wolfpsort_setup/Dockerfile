# Isolated WoLF PSORT 0.2 environment.
#
# The distribution ships NO source for its one compiled program (bin/wolfPredict);
# only prebuilt binaries for Linux/i386 and Solaris/SPARC. The host here is
# macOS/arm64, which cannot execute either. This image therefore runs the whole
# package under linux/386 emulation, which is the only route to a reproducible
# local installation of the published artefact.
#
# Nothing here touches the project's ESM/Python environment.
FROM i386/debian:bookworm-slim

# perl-base is Essential in Debian but Pod::Usage lives in perl-modules, which the
# wrappers require; install the full perl package and pin nothing else.
RUN apt-get update \
 && apt-get install -y --no-install-recommends perl \
 && rm -rf /var/lib/apt/lists/*

COPY . /opt/wolfpsort

# INSTALL step 1: copy the platform-specific binaries into ./bin.
RUN cp /opt/wolfpsort/bin/binByPlatform/binary-i386/WoLFPSORTpredictAndAlign \
       /opt/wolfpsort/bin/binByPlatform/binary-i386/runAddAlignmentToHtml \
       /opt/wolfpsort/bin/binByPlatform/binary-i386/wolfPredict \
       /opt/wolfpsort/bin/ \
 && chmod +x /opt/wolfpsort/bin/wolfPredict \
             /opt/wolfpsort/bin/WoLFPSORTpredictAndAlign \
             /opt/wolfpsort/bin/runAddAlignmentToHtml

# INSTALL step 2 (dbm feature tables; only needed by runWolfPsortHtmlTables,
# not by runWolfPsortSummary). Run it so the install is complete, and do not
# fail the build if the local dbm flavour is unavailable.
RUN cd /opt/wolfpsort/bin/psortModifiedForWolfFiles \
 && (./psortModifiedForWoLF -t all.seq > /opt/wolfpsort/install_step2.log 2>&1 \
     && echo OK >> /opt/wolfpsort/install_step2.log \
     || echo "STEP2_FAILED (runWolfPsortSummary is unaffected)" >> /opt/wolfpsort/install_step2.log) \
 && ls -la /opt/wolfpsort/bin/psortModifiedForWolfFiles >> /opt/wolfpsort/install_step2.log

# The wrappers test $ENV{LANG} against /^ja_JP/ without checking it is defined,
# which emits an uninitialised-value warning on stderr when LANG is unset.
# Setting it to C silences that and pins the message locale.
ENV LANG=C

WORKDIR /opt/wolfpsort
ENTRYPOINT ["/opt/wolfpsort/bin/runWolfPsortSummary"]

# Licensing Axiolid

Axiolid is dual-licensed: **`MPL-2.0 OR LicenseRef-Commercial`**. You may use it under the [Mozilla Public License 2.0](https://www.mozilla.org/MPL/2.0/), or under a separate commercial licence from the copyright holder ([ADR 0089](../adr/0089-dual-licence-with-a-commercial-option.md)). The canonical texts are [LICENSE](https://github.com/axiolid/kernel/blob/main/LICENSE) and [LICENSE-COMMERCIAL](https://github.com/axiolid/kernel/blob/main/LICENSE-COMMERCIAL).

One package is the exception. `axiolid-mesh-boolean-boolmesh` contains absorbed upstream code with a third-party copyright holder ([ADR 0047](../adr/0047-absorb-mesh-boolean.md)), so it is **MPL-2.0 only**.

## Under MPL-2.0

MPL is *file-level copyleft*. It keeps the kernel open while letting a product use it without relicensing the product itself.

- **Axiolid source:** If you distribute a modified Axiolid source file, make that file’s source available to recipients under MPL-2.0.
- **Your application:** Files you add for an application may remain proprietary, including when the application depends on or links Axiolid.
- **A larger work:** Distributing an executable or product that combines Axiolid with other code does not by itself place the other code under MPL-2.0.

When you ship a binary that includes Axiolid, retain the copyright and licence notices and make the corresponding source for the MPL-covered Axiolid files available to recipients. A practical approach is a `third-party-notices` page plus a source-offer URL or archive for the exact Axiolid revision and any patches.

If you modify only application files, those files are outside Axiolid’s MPL-covered source boundary. If you copy Axiolid code into an application file or modify an Axiolid file, that source-file obligation applies to the resulting file.

This avoids LGPL’s usual static-linking/relinking friction in Rust while still requiring distributed improvements to the kernel itself to stay open.

## Under a commercial licence

A commercial licence is a written agreement with the copyright holder. Its terms are those of the agreement, for example covering modified Axiolid files that you do not want to publish. `LICENSE-COMMERCIAL` only records that the option exists; it grants nothing by itself. To ask for one, open a discussion or contact the maintainer through the [repository](https://github.com/axiolid/kernel).

The commercial option never covers third-party code that Axiolid contains under its original licence, such as the absorbed mesh boolean.

## Contributor terms

The commercial option only works if the copyright holder may license every part of Axiolid that way. By submitting a contribution (a pull request, patch or other code or documentation intended for inclusion), you agree to these terms:

1. **You keep your copyright.** You remain the owner of your contribution.
2. **Licence to the project.** You license your contribution under MPL-2.0, and you also grant the Axiolid copyright holder a perpetual, worldwide, non-exclusive, royalty-free, irrevocable licence to use, reproduce, modify, distribute and sublicense it, and to license it under other terms, including the commercial licence.
3. **Patents.** You grant the recipients of Axiolid a patent licence on the same terms as section 2 of MPL-2.0, for the patent claims you own that your contribution practises.
4. **You may make the grant.** The contribution is your own work, or you have the right to submit it under these terms (for example, your employer allows it), and you name any third-party code it contains together with its licence.
5. **Open availability stays.** Every contribution accepted into Axiolid stays available under MPL-2.0.

Keep third-party code in its own dependency or clearly documented import path. Do not copy code into Axiolid without recording its provenance and licence compatibility. Third-party code is not covered by your grant, and it is not available under the commercial option unless its own licence allows that.

## FAQ

### Can a closed-source application use Axiolid?

Yes. MPL-2.0 permits this. The application’s separate source files do not become MPL-2.0 merely because they use Axiolid.

### Must kernel changes be open sourced?

Under MPL-2.0, when you distribute those modified files, yes: provide their source under MPL-2.0 to the recipients of the binary or source distribution. A commercial licence can set other terms.

### Is this legal advice?

No. This page explains the project’s intended licensing model, not legal advice. Consult counsel for a product-specific distribution decision.

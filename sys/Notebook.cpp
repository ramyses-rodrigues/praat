/* Notebook.cpp
 *
 * Copyright (C) 2023 Paul Boersma
 *
 * This code is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This code is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this work. If not, see <http://www.gnu.org/licenses/>.
 */

#include "Notebook.h"

Thing_implement (Notebook, SimpleString, 0);

autoNotebook Notebook_createFromFile (MelderFile file) {
	autoNotebook me = Thing_new (Notebook);
	my string = Melder_dup (MelderFile_peekPath (file));
	return me;
}

/* End of file Notebook.cpp */
